package cli.utils.general;

import cli.Main;
import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.tools.MatrixTools;

import java.util.ArrayList;
import java.util.List;

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


    public static List<Pixel> getAllEnrichedPixels(float[][] image, double threshold,
                                                   Zscore zscore) {
        List<Pixel> pixels = new ArrayList<>();
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                if (image[i][j] > threshold) {
                    pixels.add(new Pixel(i, j, image[i][j], (float) zscore.getZscore(image[i][j])));
                }
            }
        }
        return pixels;
    }

    public static Zscore getZscore(float[][] image, int minVal) {
        Welford welford = new Welford();
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[i].length; j++) {
                if (image[i][j] > minVal) {
                    welford.addValue(image[i][j]);
                }
            }
        }
        return welford.getZscore();
    }

    public static void saveIfVerbose(String s, float[][] array) {
        if (Main.printVerboseComments) {
            MatrixTools.saveMatrixTextNumpy(s, array);
        }

    }
}
