package cli.utils.pinpoint;

public class ArrayTools {

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
}
