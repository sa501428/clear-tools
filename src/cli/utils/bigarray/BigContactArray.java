/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2020-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package cli.utils.bigarray;

import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.datastructures.ListOfIntArrays;

import java.util.ArrayList;
import java.util.List;

public class BigContactArray {

    protected final List<int[]> binXs = new ArrayList<>();
    protected final List<int[]> binYs = new ArrayList<>();
    protected final List<float[]> binVals = new ArrayList<>();
    private final long matrixSize;

    public BigContactArray(long matrixSize) {
        this.matrixSize = matrixSize;
    }


    public void addSubList(int[] x, int[] y, float[] c) {
        binXs.add(x);
        binYs.add(y);
        binVals.add(c);
    }

    public void addSubList(int[] x, int[] y, float[] c, int counter) {
        int[] x2 = new int[counter];
        int[] y2 = new int[counter];
        float[] c2 = new float[counter];
        System.arraycopy(x, 0, x2, 0, counter);
        System.arraycopy(y, 0, y2, 0, counter);
        System.arraycopy(c, 0, c2, 0, counter);
        addSubList(x2, y2, c2);
    }

    public void addAllSubLists(BigContactArray other) {
        binXs.addAll(other.binXs);
        binYs.addAll(other.binYs);
        binVals.addAll(other.binVals);
    }

    public void clear() {
        binXs.clear();
        binYs.clear();
        binVals.clear();
    }

    private int getNumThreads() {
        return Math.min(8, binXs.size()); // todo
    }

    public long getMatrixSize() {
        return matrixSize;
    }

    public ListOfFloatArrays normalizeVectorByScaleFactor(ListOfFloatArrays newNormVector) {
        for (long k = 0; k < newNormVector.getLength(); k++) {
            float kVal = newNormVector.get(k);
            if (kVal <= 0 || Double.isNaN(kVal)) {
                newNormVector.set(k, Float.NaN);
            } else {
                newNormVector.set(k, 1.f / kVal);
            }
        }

        double normalizedSumTotal = 0, sumTotal = 0;

        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            float[] subBinVals = binVals.get(sIndx);

            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                float counts = subBinVals[z];

                double valX = newNormVector.get(x);
                double valY = newNormVector.get(y);

                if (!Double.isNaN(valX) && !Double.isNaN(valY)) {
                    double normalizedValue = counts / (valX * valY);
                    normalizedSumTotal += normalizedValue;
                    sumTotal += counts;
                    if (x != y) {
                        normalizedSumTotal += normalizedValue;
                        sumTotal += counts;
                    }
                }
            }
        }

        double scaleFactor = Math.sqrt(normalizedSumTotal / sumTotal);
        newNormVector.multiplyEverythingBy(scaleFactor);
        return newNormVector;
    }

    public ListOfIntArrays getNumNonZeroInRows() {
        ListOfIntArrays numNonZero = new ListOfIntArrays(matrixSize, 0);
        for (int sIndx = 0; sIndx < binXs.size(); sIndx++) {
            int[] subBinXs = binXs.get(sIndx);
            int[] subBinYs = binYs.get(sIndx);
            for (int z = 0; z < subBinXs.length; z++) {
                int x = subBinXs[z];
                int y = subBinYs[z];
                numNonZero.addTo(x, 1);
                if (x != y) {
                    numNonZero.addTo(y, 1);
                }
            }
        }
        return numNonZero;
    }
}
