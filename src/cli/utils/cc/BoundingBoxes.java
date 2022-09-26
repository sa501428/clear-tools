package cli.utils.cc;

import cli.Main;
import cli.utils.general.Utils;
import cli.utils.pinpoint.ArrayTools;
import cli.utils.pinpoint.ConvolutionTools;
import javastraw.reader.block.ContactRecord;
import javastraw.tools.MatrixTools;

import java.util.*;

public class BoundingBoxes {

    public static List<int[]> getEnrichedRegions(List<ContactRecord> records, int binXStart, int binYStart,
                                                 int matrixWidth, int scalar, String saveString, float[][] kernel) {
        float[][] outputD10 = new float[matrixWidth / scalar][matrixWidth / scalar];
        Utils.fillInMatrixFromRecords(outputD10, records, binXStart, binYStart, scalar);
        if (Main.printVerboseComments) {
            MatrixTools.saveMatrixTextNumpy(saveString + ".raw.S10.npy", outputD10);
        }
        float[][] kde = ConvolutionTools.sparseConvolution(outputD10, kernel);
        outputD10 = null;
        if (Main.printVerboseComments) {
            MatrixTools.saveMatrixTextNumpy(saveString + ".kde.S10.npy", kde);
        }
        float threshold = 0.8f * ArrayTools.getMax(kde);
        //Zscore zscore = LandScape.getZscore(kde, 1);
        List<Pixel> enrichedPixels = LandScape.getAllEnrichedPixels(kde, threshold); //zscore);
        return getBoundingBoxes(enrichedPixels, 2, scalar, kde.length, kde[0].length);
    }

    public static List<int[]> getBoundingBoxes(List<Pixel> enrichedPixels, int radius, int scalar,
                                               int rowLimit, int colLimit) {
        List<Pixel> sortedPixels = new ArrayList<>(enrichedPixels);
        sortedPixels.sort((o1, o2) -> Float.compare(o1.value, o2.value));
        Collections.reverse(sortedPixels);
        List<int[]> bounds = new ArrayList<>();

        while (!sortedPixels.isEmpty()) {
            Pixel pixel = sortedPixels.get(0);
            if (pixel != null) {
                int minR = pixel.row - radius;
                int minC = pixel.col - radius;
                int maxR = pixel.row + radius + 1;
                int maxC = pixel.col + radius + 1;

                Set<Pixel> toRemove = new HashSet<>();
                toRemove.add(pixel);

                while (toRemove.size() > 0) {
                    sortedPixels.removeAll(toRemove);
                    toRemove.clear();
                    for (Pixel px : sortedPixels) {
                        if (Pixel.contains(px, minR, maxR, minC, maxC)) {
                            toRemove.add(px);
                            minR = Math.min(minR, px.row - radius);
                            minC = Math.min(minC, px.col - radius);
                            maxR = Math.max(maxR, px.row + radius + 1);
                            maxC = Math.max(maxC, px.col + radius + 1);
                        }
                    }
                }

                minR = Math.max(minR, 0);
                minC = Math.max(minC, 0);
                maxR = Math.min(maxR, rowLimit);
                maxC = Math.min(maxC, colLimit);

                bounds.add(new int[]{minR * scalar, minC * scalar, maxR * scalar, maxC * scalar});
            }
        }
        return bounds;
    }
}
