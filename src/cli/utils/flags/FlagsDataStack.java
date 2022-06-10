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

import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.io.File;

public class FlagsDataStack {

    private final File dataDirectory;
    private final String customPrefix;
    private final int n;
    private final double[][] apaMatrix;

    public FlagsDataStack(int n, File outputFolder, String customPrefix) {
        this.n = n;
        apaMatrix = new double[n][n];
        dataDirectory = outputFolder;
        this.customPrefix = customPrefix;
        HiCFileTools.createValidDirectory(dataDirectory.getAbsolutePath());
    }

    public void exportData() {
        FlagsRegionStatistics stats = new FlagsRegionStatistics(apaMatrix);
        MatrixTools.saveMatrixTextNumpy((new File(dataDirectory, customPrefix + "apa.npy")).getAbsolutePath(),
                apaMatrix);
        MatrixTools.saveMatrixTextNumpy((new File(dataDirectory, customPrefix + "stats.npy")).getAbsolutePath(),
                stats.getAllValues());
    }

    public void addData(double[][] newData) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                apaMatrix[i][j] += newData[i][j];
            }
        }
    }
}
