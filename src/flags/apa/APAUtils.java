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

package flags.apa;


import javastraw.feature2D.Feature2D;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;

/**
 * Created by Muhammad Shamim on 1/21/15.
 */
public class APAUtils {

    public static RealMatrix standardNormalization(RealMatrix matrix) {
        return matrix.copy().scalarMultiply(1. / Math.max(1., mean(matrix)));
    }

    public static double mean(RealMatrix matrix) {
        return APARegionStatistics.statistics(matrix.getData()).getMean();
    }

    public static void addLocalizedData(double[][] matrix, MatrixZoomData zd, Feature2D loop,
                                        int L, int resolution, int window, NormalizationType norm, final Object key) {
        long loopX = loop.getMidPt1() / resolution;
        long loopY = loop.getMidPt2() / resolution;
        long binXStart = loopX - window;
        long binXEnd = loopX + (window + 1);
        long binYStart = loopY - window;
        long binYEnd = loopY + (window + 1);

        addLocalBoundedRegion(matrix, zd, binXStart, binXEnd, binYStart, binYEnd, L, norm, key);
    }

    public static void addLocalBoundedRegion(double[][] matrix, MatrixZoomData zd, long binXStart, long binXEnd,
                                             long binYStart, long binYEnd, int dim,
                                             NormalizationType normalizationType, final Object key) {
        List<Block> blocks;
        synchronized (key) {
            blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, normalizationType, false);
        }
        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) {
                            // only called for small regions - should not exceed int
                            int relativeX = (int) (rec.getBinX() - binXStart);
                            int relativeY = (int) (rec.getBinY() - binYStart);
                            if (relativeX >= 0 && relativeX < dim) {
                                if (relativeY >= 0 && relativeY < dim) {
                                    matrix[relativeX][relativeY] += rec.getCounts();
                                }
                            }
                        }
                    }
                }
            }
        }
        // force cleanup
        blocks.clear();
        blocks = null;
        //System.gc();
    }
}
