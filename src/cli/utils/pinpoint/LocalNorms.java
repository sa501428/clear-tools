package cli.utils.pinpoint;

public class LocalNorms {
    public static void normalizeLocally(float[][] output) {
        // todo
        boolean doSqrt = false; //norm.getLabel().toLowerCase().contains("sqrt");
        float[] rowMeans = new float[output.length];
        float[] colMeans = new float[output[0].length];

        int[] rowCounts = new int[output.length];
        int[] colCounts = new int[output[0].length];

        for (int i = 0; i < rowMeans.length; i++) {
            for (int j = 0; j < colMeans.length; j++) {
                if (output[i][j] > 0) {
                    rowMeans[i] += output[i][j];
                    colMeans[j] += output[i][j];
                    rowCounts[i]++;
                    colCounts[j]++;
                }
            }
        }

        divide(rowMeans, rowCounts);
        divide(colMeans, colCounts);

        for (int i = 0; i < rowMeans.length; i++) {
            for (int j = 0; j < colMeans.length; j++) {
                if (output[i][j] > 0) {
                    float nv = rowMeans[i] * colMeans[j];
                    if (doSqrt) nv = (float) Math.sqrt(nv);
                    output[i][j] = output[i][j] / nv;
                }
            }
        }
    }

    private static void divide(float[] means, int[] counts) {
        for (int k = 0; k < means.length; k++) {
            if (counts[k] > 0) {
                means[k] = means[k] / counts[k];
            }
        }
    }
}
