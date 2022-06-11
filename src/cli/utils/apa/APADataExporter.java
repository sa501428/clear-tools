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

package cli.utils.apa;


import javastraw.tools.MatrixTools;

import java.io.File;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by muhammadsaadshamim on 5/1/15.
 */
public class APADataExporter {

    public static void exportGenomeWideData(AtomicInteger[] gwPeakNumbers, File dataDirectory,
                                            boolean useAgNorm, float[][] globalOutput,
                                            double[] globalRowSum, double[] globalColSum) {

        Integer[] peakNumbers = {gwPeakNumbers[0].get(), gwPeakNumbers[1].get(), gwPeakNumbers[2].get()};
        String title = "N=" + peakNumbers[0] + "_(filtered)_" + peakNumbers[1] + "_(unique)_" +
                peakNumbers[2] + "_(total)";

        if (useAgNorm) {
            normalizeBySum(globalRowSum);
            normalizeBySum(globalColSum);
            float[][] gwAggNormAPAMatrix = normedCopyFloats(globalOutput, globalRowSum, globalColSum);
            MatrixTools.saveMatrixTextNumpy((new File(dataDirectory, "gw_apa_agg_norm_" + title + ".npy")).getAbsolutePath(),
                    gwAggNormAPAMatrix);
        } else {
            MatrixTools.saveMatrixTextNumpy((new File(dataDirectory, "gw_apa_" + title + ".npy")).getAbsolutePath(),
                    globalOutput);
        }
    }

    private static void normalizeBySum(double[] globalSum) {
        double average = getAverage(globalSum);
        if (average > 0) {
            for (int k = 0; k < globalSum.length; k++) {
                globalSum[k] /= average;
            }
        }
    }

    public static float[][] normedCopyFloats(float[][] original, double[] v1, double[] v2) {
        int n = original.length;
        float[][] matrix = new float[n][n];

        for (int r = 0; r < n; ++r) {
            for (int c = 0; c < n; ++c) {
                double normVal = (v1[r] * v2[c]);
                if (normVal > 0.0) {
                    matrix[r][c] = (float) (original[r][c] / normVal);
                } else {
                    matrix[r][c] = 0.0F;
                }
            }
        }

        return matrix;
    }

    public static double getAverage(double[] data) {
        double average = 0;
        if (data.length > 0) {
            double total = 0;
            for (double val : data) {
                total += val;
            }
            average = (total / data.length);
        }
        return average;
    }

}
