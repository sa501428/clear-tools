package cli.utils;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public class ExpectedUtils {
    public static double[] calculateExpected(MatrixZoomData zd, NormalizationType norm, int maxBinDist,
                                             boolean useLog) {
        double[] expected = new double[maxBinDist];
        long[] counts = new long[maxBinDist];

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            int dist = getDist(record);
            if (dist < maxBinDist) {
                if (useLog) {
                    expected[dist] += Math.log(1 + record.getCounts());
                } else {
                    expected[dist] += record.getCounts();
                }
                counts[dist]++;
            }
        }

        for (int z = 0; z < maxBinDist; z++) {
            if (counts[z] > 0) {
                expected[z] /= counts[z];
            }
        }

        if (useLog) {
            for (int z = 0; z < maxBinDist; z++) {
                expected[z] = Math.expm1(expected[z]);
            }
        }

        return expected;
    }

    public static int getDist(Feature2D loop, int resolution) {
        int binXStart = (int) (loop.getMidPt1() / resolution);
        int binYStart = (int) (loop.getMidPt2() / resolution);
        return Math.abs(binXStart - binYStart);
    }

    public static int getDist(ContactRecord record) {
        return Math.abs(record.getBinX() - record.getBinY());
    }
}
