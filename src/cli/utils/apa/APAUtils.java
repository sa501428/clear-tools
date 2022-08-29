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

package cli.utils.apa;

import javastraw.feature2D.Feature2D;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Muhammad Shamim on 1/21/15.
 */
public class APAUtils {
    public static RealMatrix standardNormalization(RealMatrix matrix) {
        return matrix.copy().scalarMultiply(1. /
                Math.max(1., statistics(matrix.getData()).getMean()));
    }

    public static DescriptiveStatistics statistics(double[][] x) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double[] row : x)
            for (double val : row)
                stats.addValue(val);
        return stats;
    }


    public static RealMatrix centerNormalization(RealMatrix matrix) {

        int center = matrix.getRowDimension() / 2;
        double centerVal = matrix.getEntry(center, center);

        if (centerVal == 0) {
            centerVal = MatrixTools.minimumPositive(matrix.getData());
            if (centerVal == 0)
                centerVal = 1;
        }

        return matrix.copy().scalarMultiply(1. / centerVal);
    }

    public static double peakEnhancement(RealMatrix matrix) {
        int rows = matrix.getRowDimension();
        int center = rows / 2;
        double centerVal = matrix.getEntry(center, center);
        double remainingSum = MatrixTools.sum(matrix.getData()) - centerVal;
        double remainingAverage = remainingSum / (rows * rows - 1);
        return centerVal / remainingAverage;
    }

    public static RealMatrix rankPercentile(RealMatrix data) {
        int n = data.getColumnDimension();
        StatPercentile percentile = new StatPercentile(MatrixTools.flattenedRowMajorOrderMatrix(data.getData()));

        RealMatrix matrix = new Array2DRowRealMatrix(n, n);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                double currValue = data.getEntry(r, c);
                if (currValue == 0) {
                    matrix.setEntry(r, c, 0);
                } else {
                    matrix.setEntry(r, c, percentile.evaluate(currValue));
                }
                //matrix.setEntry(r, c, percentile.evaluate());
            }
        }
        return matrix;
    }

    public static ArrayList<Feature2D> filterFeaturesBySize(List<Feature2D> features,
                                                            double minPeakDist, double maxPeakDist, int resolution) {
        ArrayList<Feature2D> sizeFilteredFeatures = new ArrayList<>();

        for (Feature2D feature : features) {
            double xMidPt = feature.getMidPt1();
            double yMidPt = feature.getMidPt2();
            int dist = (int) Math.round(Math.abs(xMidPt - yMidPt) / resolution);

            if (dist >= minPeakDist) {
                if (dist <= maxPeakDist) {
                    sizeFilteredFeatures.add(feature);
                }
            }
        }
        return new ArrayList<>(sizeFilteredFeatures);
    }

    public static void inPlaceSumVectors(double[] globalSum, double[] vector) {
        for (int j = 0; j < globalSum.length; j++) {
            if (vector[j] > 0) {
                globalSum[j] += vector[j];
            }
        }
    }

    public static void inPlaceSumMatrices(float[][] globalSum, float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > 0) {
                    globalSum[i][j] += matrix[i][j];
                }
            }
        }
    }

    public static void addLocalRowSums(double[] sums, double[] vector, int binStart) {
        for (int i = 0; i < sums.length; i++) {
            double val = vector[binStart + i];
            if (val > 0) {
                sums[i] += val;
            }
        }
    }
}
