package cli.utils.clique;

import javastraw.tools.ParallelizationTools;

import java.util.concurrent.atomic.AtomicInteger;

public class MatrixUtils {

    public static float[][] cube(float[][] a) {
        return multiply(a, multiply(a, a));
    }

    public static float[][] multiply(float[][] a, float[][] b) {
        int numARows = a.length;
        int numACols = a[0].length;
        int numBCols = b[0].length;

        float[][] result = new float[numARows][numBCols];
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < numARows) {
                for (int j = 0; j < numBCols; j++) {
                    for (int k = 0; k < numACols; k++) {
                        result[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        });
        return result;
    }

    public static void setEverythingBelowDiagonalToZero(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j <= i; j++) {
                matrix[i][j] = 0;
            }
        }
    }
}
