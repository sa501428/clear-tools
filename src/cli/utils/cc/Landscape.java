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

/**
 * Implementation of 2-pass algorithm for finding connected components
 */
public class Landscape {
/*
    private static final int NOT_SET = 0;
    private static final int IN_QUEUE = -1;

    public static final float ABSOLUTE_CUTOFF = 10f;

    public static void extractMaxima(float[][] kde, int binXStart, int binYStart, long resolution,
                                     List<Feature2D> pinpointedLoops, Feature2D loop, String saveString,
                                     boolean onlyGetOne) {

        Pixel max = getMax(kde, (int) Math.max(100 / resolution, 3));
        if (max != null) {
            long start1 = resolution * (binXStart + max.x);
            long start2 = resolution * (binYStart + max.y);
            long end1 = start1 + resolution;
            long end2 = start2 + resolution;

            Map<String, String> map = loop.getAttributes();
            map.put("threshold", "" + kde[max.x][max.y]);

            Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                    loop.getChr2(), start2, end2, Color.BLACK, map);
            pinpointedLoops.add(feature);
        }

        {
            float threshold = ArrayTools.getMax(kde) * 0.85f;
            if (threshold > ABSOLUTE_CUTOFF) {
                List<LocalMaxima> maxima = detect(kde, threshold, saveString);

                if (maxima.size() > 0) {
                    LocalMaxima best = maxima.get(0);

                    for (LocalMaxima max : maxima) {
                        Map<String, String> attributes = new HashMap<>();
                        attributes.put("pinpoint_area", "" + max.area);
                        attributes.put("pinpoint_max_enrichment", "" + max.maxVal);
                        attributes.put("pinpoint_mean_enrichment", "" + max.meanVal);

                        if (best.maxVal < max.maxVal) {
                            best = max;
                        }

                        long start1 = resolution * (binXStart + max.maxCoordinate.x);
                        long start2 = resolution * (binYStart + max.maxCoordinate.y);
                        long end1 = start1 + resolution;
                        long end2 = start2 + resolution;

                        Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                                loop.getChr2(), start2, end2, Color.BLACK, attributes);
                        pinpointedLoops.add(feature);
                    }

                    best = null;
                    maxima.clear();
                }
            }
        }
    }

    private static List<Feature2D> removeOverlappingPixels(List<Feature2D> features, boolean useNMS) {
        Map<SimpleLocation, LinkedList<Feature2D>> map = groupNearbyRecords(simpleFilter(features), 1000000);
        Set<Feature2D> coalesced = new HashSet<>();

        for (LinkedList<Feature2D> featureLL : map.values()) {
            featureLL.sort((o1, o2) -> {
                int val = Long.compare(o1.getWidth1(), o2.getWidth1());
                if (val == 0) {
                    return Long.compare(o1.getWidth2(), o2.getWidth2());
                }
                return val;
            });

            while (!featureLL.isEmpty()) {

                Feature2D pixel = featureLL.pollFirst();
                if (pixel != null) {
                    featureLL.remove(pixel);
                    int buffer = (int) Math.max(pixel.getWidth1(), pixel.getWidth2());

                    Set<Feature2D> pixelList = new HashSet<>();
                    pixelList.add(pixel);
                    for (Feature2D px : featureLL) {
                        if (hasOverlap(px, pixel, buffer)) {
                            pixelList.add(px);
                        }
                    }

                    featureLL.removeAll(pixelList);
                    pixel.addStringAttribute("NumCollapsed", String.valueOf(pixelList.size()));

                    long start1 = FeatureStats.minStart1(pixelList);
                    long start2 = FeatureStats.minStart2(pixelList);
                    long end1 = FeatureStats.maxEnd1(pixelList);
                    long end2 = FeatureStats.maxEnd2(pixelList);
                    coalesced.add(new Feature2D(Feature2D.FeatureType.PEAK, pixel.getChr1(), start1, end1,
                            pixel.getChr2(), start2, end2, Color.BLACK, pixel.getAttributes()));

                    pixelList.clear();
                }
            }
        }
        map.clear();

        return new ArrayList<>(coalesced);
    }

    public static List<LocalMaxima> detect(float[][] image, double threshold, String saveString) {
        int r = image.length;
        int c = image[0].length;

        int[][] labels = new int[r][c];
        int nextLabel = 1;

        List<LocalMaxima> results = new ArrayList<>();
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                if (image[i][j] > threshold && labels[i][j] == NOT_SET) {
                    Queue<Pixel> points = new LinkedList<>();
                    points.add(new Pixel(i, j));
                    labels[i][j] = IN_QUEUE; // actively being processed
                    LocalMaxima maxima = processRegion(image, threshold, labels, points, nextLabel);
                    results.add(maxima);
                    nextLabel++;
                }
            }
        }

        labels = null;
        return results;
    }

    public static Pixel getMax(float[][] image, int dist) {
        double maxVal = ArrayTools.getMax(image);
        if (maxVal > 0) {
            LinkedList<Pixel> pixels = getAllEnrichedPixels(image, maxVal * 0.99);
            if (pixels.size() == 1) {
                return pixels.get(0);
            }
            return pickBestCentroid(pixels, dist);
        }
        return null;
    }

    public static Pixel pickBestCentroid(LinkedList<Pixel> regions, int dist) {
        List<PixelSet> pixelSets = new ArrayList<>();
        while (regions.size() > 0) {
            Pixel pixel = regions.pollFirst();
            PixelSet pixelSet = new PixelSet(pixel);
            regions.remove(pixel);

            int prevSize = 0;
            while (prevSize != pixelSet.size()) {
                prevSize = pixelSet.size();
                for (Pixel px : regions) {
                    if (pixelSet.distanceTo(px) < dist) {
                        pixelSet.add(px);
                    }
                }
                regions.removeAll(pixelSet.pixels);
            }
            if (pixelSet.size() > 1) { // centroids should have more than 1
                pixelSets.add(pixelSet);
            }
        }

        if (pixelSets.size() > 0) {
            pixelSets.sort(Comparator.comparingInt(PixelSet::size));
            Collections.reverse(pixelSets);
            return pixelSets.get(0).toPixel();
        }
        return null;
    }

    private static LinkedList<Pixel> getAllEnrichedPixels(float[][] image, double threshold) {
        LinkedList<Pixel> pixels = new LinkedList<>();
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                if (image[i][j] > threshold) {
                    pixels.add(new Pixel(i, j));
                }
            }
        }
        return pixels;
    }

    private static LocalMaxima processRegion(float[][] image, double threshold, int[][] labels, Queue<Pixel> points, int id) {
        if (points.isEmpty()) return null;

        int area = 0;
        Pixel maxCoordinate = points.peek();
        float maxVal = image[maxCoordinate.x][maxCoordinate.y];

        while (points.size() > 0) {
            Pixel current = points.poll();
            if (current != null) {
                area++;
                labels[current.x][current.y] = id;

                if (image[current.x][current.y] > maxVal) {
                    maxVal = image[current.x][current.y];
                    maxCoordinate = current;
                }

                for (int i = current.x - 1; i < current.x + 2; i++) {
                    if (i < 0 || i >= image.length) continue;
                    for (int j = current.y - 1; j < current.y + 2; j++) {
                        if (j < 0 || j >= image[i].length) continue;
                        if (image[i][j] > threshold && labels[i][j] == NOT_SET) {
                            points.add(new Pixel(i, j));
                            labels[i][j] = IN_QUEUE; // actively being processed
                        }
                    }
                }
            }
        }

        return new LocalMaxima(maxCoordinate, maxVal, area);
    }
    
 */
}
