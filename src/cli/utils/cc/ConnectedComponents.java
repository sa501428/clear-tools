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

import cli.utils.ArrayTools;
import javastraw.feature2D.Feature2D;

import java.awt.*;
import java.util.List;
import java.util.Queue;
import java.util.*;

/**
 * Implementation of 2-pass algorithm for finding connected components
 */
public class ConnectedComponents {

    private static final int NOT_SET = 0;
    private static final int IN_QUEUE = -1;

    public static void extractMaxima(float[][] kde, int binXStart, int binYStart, long resolution,
                                     List<Feature2D> pinpointedLoops, Feature2D loop, String saveString) {
        float threshold = ArrayTools.getMax(kde) * 0.85f;
        if (threshold > 10) {
            List<LocalMaxima> maxima = detect(kde, threshold, saveString);
            for (LocalMaxima max : maxima) {
                Map<String, String> attributes = new HashMap<>();
                attributes.put("pinpoint_area", "" + max.area);
                attributes.put("pinpoint_enrichment", "" + max.maxVal);

                long start1 = resolution * (binXStart + max.maxCoordinate.x);
                long start2 = resolution * (binYStart + max.maxCoordinate.y);
                long end1 = start1 + resolution;
                long end2 = start2 + resolution;

                Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, loop.getChr1(), start1, end1,
                        loop.getChr2(), start2, end2, Color.BLACK, attributes);
                pinpointedLoops.add(feature);
            }
            maxima.clear();
        }
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

        //if (Main.printVerboseComments) {
        //    MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_kde.npy")).getAbsolutePath(), image);
        //    MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_label.npy")).getAbsolutePath(), labels);
        //}

        labels = null;
        return results;
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
}
