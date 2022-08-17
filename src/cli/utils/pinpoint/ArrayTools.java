package cli.utils.pinpoint;

public class ArrayTools {

    public static float getSum(float[][] matrix) {
        float total = 0;
        for (float[] row : matrix) {
            for (float val : row) {
                total += val;
            }
        }
        return total;
    }

    public static float getMax(float[][] matrix) {
        float max = matrix[0][0];
        for (float[] row : matrix) {
            for (float val : row) {
                if (val > max) {
                    max = val;
                }
            }
        }
        return max;
    }

    public static void inPlaceDivideArrayBy(float[][] result, float val) {
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[i].length; j++) {
                result[i][j] /= val;
            }
        }
    }
}
