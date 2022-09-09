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


import javastraw.expected.ExpectedModel;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.List;

/**
 * Created by Muhammad Shamim on 1/21/15.
 */
public class Utils {

    public static void addLocalizedData(float[][] matrix, MatrixZoomData zd, Feature2D loop,
                                        int matrixWidth, int resolution, int window, NormalizationType norm) {
        long binXStart = (loop.getMidPt1() / resolution) - window;
        long binYStart = (loop.getMidPt2() / resolution) - window;
        addLocalBoundedRegion(matrix, zd, binXStart, binYStart, matrixWidth, norm);
    }

    public static void addLocalBoundedRegion(float[][] matrix, MatrixZoomData zd, long binXStart, long binYStart,
                                             int matrixWidth, NormalizationType norm) {

        long binXEnd = binXStart + (matrixWidth + 1);
        long binYEnd = binYStart + (matrixWidth + 1);
        List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart,
                binXEnd, binYEnd, norm, false);

        fillInMatrixFromBlocks(matrix, blocks, binXStart, binYStart);
        blocks.clear();
        blocks = null;
    }

    public static float[][] getRegion(MatrixZoomData zd, long binXStart, long binYStart,
                                      long binXEnd, long binYEnd, NormalizationType norm) {
        int numRows = (int) (binXEnd - binXStart);
        int numCols = (int) (binYEnd - binYStart);
        float[][] matrix = new float[numRows][numCols];
        List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart,
                binXEnd, binYEnd, norm, false);
        fillInMatrixFromBlocks(matrix, blocks, binXStart, binYStart);
        blocks.clear();
        blocks = null;
        return matrix;
    }

    public static void fillInMatrixFromBlocks(float[][] matrix, List<Block> blocks, long binXStart, long binYStart) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;

        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) {
                            // only called for small regions - should not exceed int
                            int relativeX = (int) (rec.getBinX() - binXStart);
                            int relativeY = (int) (rec.getBinY() - binYStart);
                            if (relativeX >= 0 && relativeX < numRows) {
                                if (relativeY >= 0 && relativeY < numCols) {
                                    matrix[relativeX][relativeY] += rec.getCounts();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public static void fillInOEMatrixFromBlocks(float[][] matrix, List<Block> blocks,
                                                long binXStart, long binYStart,
                                                ExpectedValueFunction df, int chrIndex, double pseudocount) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;

        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) {
                            // only called for small regions - should not exceed int

                            int relativeX = (int) (rec.getBinX() - binXStart);
                            int relativeY = (int) (rec.getBinY() - binYStart);
                            if (relativeX >= 0 && relativeX < numRows) {
                                if (relativeY >= 0 && relativeY < numCols) {
                                    int dist = Math.abs(rec.getBinX() - rec.getBinY());
                                    double expected = df.getExpectedValue(chrIndex, dist);
                                    matrix[relativeX][relativeY] = (float) ((rec.getCounts() + pseudocount) /
                                            (expected + pseudocount));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public static void fillInExpectedMatrix(float[][] matrix, Feature2D loop,
                                            int matrixWidth, ExpectedModel expectedVector,
                                            int resolution, int window) {

        long binXStart = (loop.getMidPt1() / resolution) - window;
        long binYStart = (loop.getMidPt2() / resolution) - window;

        for (int relativeX = 0; relativeX < matrixWidth; relativeX++) {
            for (int relativeY = 0; relativeY < matrixWidth; relativeY++) {
                long X = relativeX + binXStart;
                long Y = relativeY + binYStart;
                int dist = (int) Math.abs(X - Y);
                double expected = expectedVector.getExpectedFromUncompressedBin(dist);
                matrix[relativeX][relativeY] = (float) expected;
            }
        }
    }
}
