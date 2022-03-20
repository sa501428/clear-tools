package cli.apa;

public class AggNorm {

    public static double[][] getAggNormedMatrix(double[][] original){
        double[] rowSums = new double[original.length];
        double[] colSums = new double[original[0].length];
        for(int r = 0; r < original.length; ++r) {
            for(int c = 0; c < original[0].length; ++c) {
                if (original[r][c] > 0) {
                    rowSums[r] += original[r][c];
                    colSums[c] += original[r][c];
                }
            }
        }
        scalarMultiply(rowSums, 1.0 / getAverage(rowSums));
        scalarMultiply(colSums, 1.0 / getAverage(colSums));
        return normedCopy(original, rowSums, colSums);
    }

    private static void scalarMultiply(double[] vector, double val) {
        for(int i = 0; i < vector.length; i++){
            vector[i] *= val;
        }
    }

    private static double getAverage(double[] vector) {
        double sum = 0;
        for(double val : vector){
            sum += val;
        }
        return sum / vector.length;
    }

    public static double[][] normedCopy(double[][] original, double[] rowSums, double[] colSums) {
        double[][] matrix = new double[original.length][original[0].length];

        for(int r = 0; r < original.length; ++r) {
            for(int c = 0; c < original[0].length; ++c) {
                double normVal = rowSums[r] * colSums[c];
                if (normVal > 0.0D) {
                    matrix[r][c] = original[r][c] / normVal;
                }
            }
        }

        return matrix;
    }
}
