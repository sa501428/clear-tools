package cli.utils.pinpoint;

import cli.utils.general.ArrayTools;

public class ConvolutionTools {

    public static float[][] sparseConvolution(float[][] image, float[][] kernel) {

        float[][] result = new float[image.length][image[0].length];

        int halfWidth = kernel.length / 2;
        int maxR = image.length;
        int maxC = image[0].length;

        for (int i = 0; i < maxR; i++) {
            for (int j = 0; j < maxC; j++) {
                if (image[i][j] > 0) {
                    for (int ki = 0; ki < kernel.length; ki++) {
                        int r = i - halfWidth + ki;
                        if (r >= 0 && r < maxR) {
                            for (int kj = 0; kj < kernel[ki].length; kj++) {
                                int c = j - halfWidth + kj;
                                if (c >= 0 && c < maxC) {
                                    result[r][c] += image[i][j] * kernel[ki][kj];
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    }

    public static float[][] getManhattanKernel(int width) {
        float[][] kernel = new float[width][width];
        int halfWidth = width / 2;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                int d2 = Math.abs(i - halfWidth) + Math.abs(j - halfWidth);
                if (d2 < halfWidth) {
                    kernel[i][j] = halfWidth - d2;
                }
            }
        }

        float maxK = ArrayTools.getMax(kernel);

        // normalize center to 1
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                kernel[i][j] /= maxK;
            }
        }
        return kernel;
    }
}
