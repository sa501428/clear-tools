package cli.utils;

public class ArrayTools {

    public static float getSum(float[][] matrix) {
        double total = 0;
        for (float[] row : matrix) {
            for (float val : row) {
                total += val;
            }
        }
        total /= (matrix.length * matrix[0].length);
        return (float) total;
    }

    public static int getSum(int[] histogram) {
        int total = 0;
        for (int val : histogram) {
            total += val;
        }
        return total;
    }

}
