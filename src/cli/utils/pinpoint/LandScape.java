/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package cli.utils.pinpoint;

import cli.utils.general.ArrayTools;
import cli.utils.general.Pixel;
import cli.utils.general.Utils;
import javastraw.expected.Zscore;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.awt.*;
import java.util.List;
import java.util.*;

public class LandScape {

    private static final int MIN_ENRICHED_PIXELS = 4;
    private static final int PEAK_WIDTH_LIMIT_1D = 10;

    public static void extractMaxima(List<ContactRecord> records, int originalBinXStart, int originalBinYStart, long resolution,
                                     List<Feature2D> pinpointedLoops, List<Feature2D> pinpointedBounds,
                                     Feature2D loop, String saveString,
                                     boolean onlyGetOne, int matrixWidth, float[][] kernel, float[][] compressedKernel) {

        List<int[]> bounds = BoundingBoxes.getEnrichedRegions(records, originalBinXStart, originalBinYStart, matrixWidth, 20,
                saveString, compressedKernel);

        int counter = 0;
        for (int[] bound : bounds) {
            if (onlyGetOne && counter > 0) return;
            // minR*scalar, minC*scalar, maxR*scalar, maxC*scalar
            int newBinXStart = originalBinXStart + bound[0];
            int newBinYStart = originalBinYStart + bound[1];
            int newMatrixWidth = Math.max(bound[2] - bound[0], bound[3] - bound[1]);

            float[][] output = new float[newMatrixWidth][newMatrixWidth];
            Utils.fillInMatrixFromRecords(output, records, newBinXStart, newBinYStart);
            ArrayTools.saveIfVerbose(saveString + ".raw.c" + counter + ".npy", output);
            float[] outputR = ArrayTools.getNormedRowSums(output);
            float[] outputC = ArrayTools.getNormedColSums(output);

            float[][] kde = ConvolutionTools.sparseConvolution(output, kernel);
            ArrayTools.saveIfVerbose(saveString + ".kde.c" + counter + ".npy", kde);
            float[] kdeR = ArrayTools.getNormedRowSums(kde);
            float[] kdeC = ArrayTools.getNormedColSums(kde);
            output = null; // clear output

            float[] rowSignal = ArrayTools.multiply(outputR, kdeR);
            float[] colSignal = ArrayTools.multiply(outputC, kdeC);

            Zscore zscore = ArrayTools.getZscore(kde, 1);
            double threshold = zscore.getValForZscore(3);
            List<Pixel> preNormEnrichments = ArrayTools.getAllEnrichedPixels(kde, threshold, zscore);

            if (preNormEnrichments.size() > MIN_ENRICHED_PIXELS) {
                LocalMaxima maxima = coalesceAndRetainMaximum(preNormEnrichments, 1,
                        rowSignal, colSignal);
                preNormEnrichments.clear();
                if (maxima != null) {
                    pinpointedBounds.add(featureFromBounds(maxima, resolution, newBinXStart, newBinYStart, loop));
                    pinpointedLoops.add(peakFromMaxima(maxima, resolution, newBinXStart, newBinYStart, loop));
                }
            }
            counter++;
        }
    }

    public static LocalMaxima coalesceAndRetainMaximum(List<Pixel> pixels, int radius,
                                                       float[] rowSignal, float[] colSignal) {
        Pixel pixel = Pixel.getMax(pixels, rowSignal, colSignal, PEAK_WIDTH_LIMIT_1D);
        if (pixel != null) {
            int numCollapsed = 0;
            int minR = pixel.row - radius;
            int minC = pixel.col - radius;
            int maxR = pixel.row + radius + 1;
            int maxC = pixel.col + radius + 1;

            Set<Pixel> toRemove = new HashSet<>();
            toRemove.add(pixel);

            while (toRemove.size() > 0) {
                pixels.removeAll(toRemove);
                numCollapsed += toRemove.size();
                toRemove.clear();

                for (Pixel px : pixels) {
                    if (Pixel.contains(px, minR, maxR, minC, maxC)) {
                        toRemove.add(px);

                        minR = Math.min(minR, px.row - radius);
                        minC = Math.min(minC, px.col - radius);
                        maxR = Math.max(maxR, px.row + radius + 1);
                        maxC = Math.max(maxC, px.col + radius + 1);
                    }
                }
            }

            if (numCollapsed > MIN_ENRICHED_PIXELS) {
                return new LocalMaxima(pixel, numCollapsed, minR, minC, maxR, maxC);
            }
        }
        return null;
    }


    private static Feature2D peakFromMaxima(LocalMaxima maxima, long resolution,
                                            int binXStart, int binYStart, Feature2D loop) {
        Pixel max = maxima.maxCoordinate;
        Map<String, String> attributes = new HashMap<>(loop.getAttributes());
        attributes.put("pinpoint_zscore", "" + max.zScore);
        attributes.put("pinpoint_value", "" + max.value);
        attributes.put("pinpoint_area", "" + maxima.area);
        attributes.put("pinpoint_start_1", "" + resolution * (binXStart + maxima.minR));
        attributes.put("pinpoint_start_2", "" + resolution * (binYStart + maxima.minC));
        attributes.put("pinpoint_end_1", "" + resolution * (binXStart + maxima.maxR));
        attributes.put("pinpoint_end_2", "" + resolution * (binYStart + maxima.maxC));

        long start1 = resolution * (binXStart + max.row);
        long start2 = resolution * (binYStart + max.col);
        long end1 = start1 + resolution;
        long end2 = start2 + resolution;

        return new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                loop.getChr2(), start2, end2, Color.BLACK, attributes);
    }

    private static Feature2D featureFromBounds(LocalMaxima maxima, long resolution,
                                               int binXStart, int binYStart, Feature2D loop) {
        Pixel max = maxima.maxCoordinate;
        Map<String, String> attributes = new HashMap<>(loop.getAttributes());
        attributes.put("pinpoint_zscore", "" + max.zScore);
        attributes.put("pinpoint_value", "" + max.value);
        attributes.put("pinpoint_area", "" + maxima.area);
        attributes.put("pinpoint_center_x", "" + resolution * (binXStart + max.row));
        attributes.put("pinpoint_center_y", "" + resolution * (binYStart + max.col));

        long start1 = resolution * (binXStart + maxima.minR);
        long start2 = resolution * (binYStart + maxima.minC);
        long end1 = resolution * (binXStart + maxima.maxR);
        long end2 = resolution * (binYStart + maxima.maxC);

        return new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                loop.getChr2(), start2, end2, Color.BLACK, attributes);
    }

}
