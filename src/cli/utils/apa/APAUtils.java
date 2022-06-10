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
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.IOException;
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

    public static RealMatrix extractLocalizedData(MatrixZoomData zd, Feature2D loop,
                                                  int L, int resolution, int window, NormalizationType norm) throws IOException {
        long loopX = loop.getMidPt1() / resolution;
        long loopY = loop.getMidPt2() / resolution;
        long binXStart = loopX - window;
        long binXEnd = loopX + (window + 1);
        long binYStart = loopY - window;
        long binYEnd = loopY + (window + 1);

        return HiCFileTools.extractLocalBoundedRegion(zd, binXStart, binXEnd, binYStart, binYEnd, L, L, norm, false);
    }

    public static List<RealMatrix> extractLocalizedRowSums(MatrixZoomData zd, Feature2D loop,
                                                           int L, int resolution, int window, NormalizationType norm) throws IOException {
        long loopX = loop.getMidPt1() / resolution;
        long loopY = loop.getMidPt2() / resolution;
        long binXStart = loopX - window;
        long binXEnd = loopX + (window + 1);
        long binYStart = loopY - window;
        long binYEnd = loopY + (window + 1);
        long chrXend = zd.getChr2().getLength() / zd.getBinSize() + 1;
        long chrYend = zd.getChr1().getLength() / zd.getBinSize() + 1;

        List<RealMatrix> vectors = new ArrayList<>();

        vectors.add(HiCFileTools.extractLocalRowSums(zd, binXStart, binXEnd, 0, chrXend, L, norm, false));
        vectors.add(HiCFileTools.extractLocalRowSums(zd, binYStart, binYEnd, 0, chrYend, L, norm, false));

        return vectors;
    }
}
