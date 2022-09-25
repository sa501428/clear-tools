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

package cli.utils.cc;

import cli.utils.pinpoint.LocalNorms;
import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.feature2D.Feature2D;

import java.awt.*;
import java.util.List;
import java.util.*;

public class LandScape {
    public static void extractMaxima(float[][] kde, int binXStart, int binYStart, long resolution,
                                     List<Feature2D> pinpointedLoops, Feature2D loop, String saveString,
                                     boolean onlyGetOne) {

        List<Pixel> enrichedPixels = getAllEnrichedPixels(kde);
        if (enrichedPixels.size() == 0) return;

        LocalNorms.normalizeLocally(kde);
        List<Pixel> normEnrichedPixels = getAllEnrichedPixels(kde);
        if (normEnrichedPixels.size() == 0) return;

        List<Pixel> maxima = twoPassCoalesceAndRetainMaxima(normEnrichedPixels,
                enrichedPixels, (int) (100 / resolution) + 1);

        enrichedPixels.clear();
        normEnrichedPixels.clear();

        for (int i = 0; i < maxima.size(); i++) {
            if (onlyGetOne && i > 0) return;

            Pixel max = maxima.get(i);
            Map<String, String> attributes = new HashMap<>(loop.getAttributes());
            attributes.put("pinpoint_zscore", "" + max.zScore);
            attributes.put("pinpoint_value", "" + max.value);

            long start1 = resolution * (binXStart + max.row);
            long start2 = resolution * (binYStart + max.col);
            long end1 = start1 + resolution;
            long end2 = start2 + resolution;

            Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                    loop.getChr2(), start2, end2, Color.BLACK, attributes);
            pinpointedLoops.add(feature);
        }
        maxima.clear();
    }

    public static List<Pixel> twoPassCoalesceAndRetainMaxima(List<Pixel> pixels, List<Pixel> enrichments, int radius) {
        List<Pixel> sortedPixels = new ArrayList<>(pixels);
        sortedPixels.sort((o1, o2) -> Float.compare(o1.value, o2.value));
        Collections.reverse(sortedPixels);

        List<Pixel> coalesced = new ArrayList<>();
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

                for (Pixel px : enrichments) {
                    if (Pixel.contains(px, minR, maxR, minC, maxC)) {
                        toRemove.add(px);
                    }
                }

                if (toRemove.size() > 0) {
                    // has enrichments in the vicinity of the region
                    enrichments.removeAll(toRemove);
                    coalesced.add(pixel);
                }
            }
        }
        return coalesced;
    }

    private static List<Pixel> getAllEnrichedPixels(float[][] image) {
        Zscore zscore = getZscore(image, 1);
        double threshold = zscore.getValForZscore(3);
        List<Pixel> pixels = new ArrayList<>();
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                if (image[i][j] > threshold) {
                    pixels.add(new Pixel(i, j, image[i][j], (float) zscore.getZscore(image[i][j])));
                }
            }
        }
        return pixels;
    }

    private static Zscore getZscore(float[][] image, int minVal) {
        Welford welford = new Welford();
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                if (image[i][j] > minVal) {
                    welford.addValue(image[i][j]);
                }
            }
        }
        return welford.getZscore();
    }
}
