package cli.utils.pinpoint;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.awt.*;
import java.util.List;
import java.util.Map;

public class CentroidFinder {

    public static void findCentroid(List<ContactRecord> records, int binXStart, int binYStart,
                                    int binXEnd, int binYEnd, long resolution,
                                    List<Feature2D> pinpointedLoops, Feature2D loop, String saveString,
                                    boolean onlyGetOne) {
        float[] cXYOld = new float[]{-1, -1};
        float[] cXY = getCentroid(records, cXYOld, 2 * Math.max(binXEnd - binXStart, binYEnd - binYStart),
                binXStart, binYStart);

        for (int radius : new int[]{500, 200}) {
            int iter = 0;
            while (distance(cXY, cXYOld) > 1e-2) {
                cXYOld = cXY;
                cXY = getCentroid(records, cXYOld, radius, binXStart, binYStart);
                if (iter++ > 100) {
                    System.err.println("Something wrong?");
                }
            }
        }

        if (cXY[0] > -1 && cXY[1] > -1) {
            long start1 = (long) (resolution * (binXStart + cXY[0]));
            long start2 = (long) (resolution * (binYStart + cXY[1]));
            long end1 = start1 + resolution;
            long end2 = start2 + resolution;

            Map<String, String> map = loop.getAttributes();

            Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                    loop.getChr2(), start2, end2, Color.BLACK, map);
            pinpointedLoops.add(feature);
        }
    }

    private static float[] getCentroid(List<ContactRecord> records, float[] previous, int radius,
                                       int binXStart, int binYStart) {
        float counts = 0;
        float totalX = 0, totalY = 0;

        for (ContactRecord record : records) {
            if (distance(previous, record, binXStart, binYStart) < radius) {
                totalX += ((record.getBinX() - binXStart) * record.getCounts());
                totalY += ((record.getBinY() - binYStart) * record.getCounts());
                counts += record.getCounts();
            }
        }

        if (counts > 0) {
            return new float[]{totalX / counts, totalY / counts};
        }
        return new float[]{-1, -1};
    }

    private static float distance(float[] a, ContactRecord b, int binXStart, int binYStart) {
        return distance(a, new float[]{b.getBinX() - binXStart, b.getBinY() - binYStart});
    }

    private static float distance(float[] a, float[] b) {
        return Math.max(Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
    }

}
