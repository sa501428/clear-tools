/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package cli.utils.flags;

import javastraw.tools.MatrixTools;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Created by muhammadsaadshamim on 5/4/15.
 */
public class FlagsRegionStatistics {
    private static final int regionWidth = 10;
    private final double peak2mean;
    private final double peak2UL;
    private final double avgUR;
    private final double peak2UR;
    private final double peak2LL;
    private final double peak2LR;
    private final double ZscoreLL;

    public FlagsRegionStatistics(double[][] data) {

        int dimension = data.length;
        int midPoint = dimension / 2;
        double centralVal = data[midPoint][midPoint];

        int dimMinusWidth = dimension - regionWidth;

        double mean = (MatrixTools.sum(data) - centralVal) / (dimension * dimension - 1);
        peak2mean = centralVal / mean;

        double avgUL = mean(getSubMatrix(data, 0, regionWidth - 1, 0, regionWidth - 1));
        peak2UL = centralVal / avgUL;

        avgUR = mean(getSubMatrix(data, 0, regionWidth - 1, dimMinusWidth, dimension - 1));
        peak2UR = centralVal / avgUR;

        double avgLL = mean(getSubMatrix(data, dimMinusWidth, dimension - 1, 0, regionWidth - 1));
        peak2LL = centralVal / avgLL;

        double avgLR = mean(getSubMatrix(data, dimMinusWidth, dimension - 1, dimMinusWidth, dimension - 1));
        peak2LR = centralVal / avgLR;

        DescriptiveStatistics yStats = statistics(getSubMatrix(data, dimMinusWidth, dimension - 1, 0, regionWidth - 1));
        ZscoreLL = (centralVal - yStats.getMean()) / yStats.getStandardDeviation();
    }

    private double[][] getSubMatrix(double[][] data, int i, int i1, int i2, int i3) {
        return new Array2DRowRealMatrix(data).getSubMatrix(i, i1, i2, i3).getData();
    }

    public static DescriptiveStatistics statistics(double[][] x) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double[] row : x)
            for (double val : row)
                stats.addValue(val);
        return stats;
    }

    private static double mean(double[][] x) {
        return statistics(x).getMean();
    }

    public double[] getAllValues() {
        return new double[]{peak2mean, avgUR, peak2UL, peak2UR, peak2LL, peak2LR, ZscoreLL};
    }
}
