package cli.utils.anchors;

public class Convolution1DTools {
    // hardcoded weights for a 3-point gaussian convolution
    // .24 .52 .24
    public static float[] smooth(float[] data) {
        float[] smooth = new float[data.length];
        smooth[0] = data[0];
        smooth[data.length - 1] = data[data.length - 1];

        for (int i = 1; i < data.length - 1; i++) {
            smooth[i] = (float) (.24 * data[i - 1] + .52 * data[i] + .24 * data[i + 1]);
        }
        return smooth;
    }

    private static final float[] weights5 = new float[]{.055f, .245f, .4f, .245f, .055f};
    // f(−2)≈0.054,f(−1)≈0.242,f(0)≈0.399,f(1)≈0.242,f(2)≈0.054
    // 0.055+0.245+0.4+0.245+0.055

    public static float[] smooth5(int[] data) {
        return smooth(data, weights5);
    }

    public static float[] smooth(int[] data, float[] weights) {
        float[] smooth = new float[data.length];

        int width = weights.length / 2;

        // copy edges
        for (int i = 0; i < width; i++) {
            smooth[i] = data[i];
        }
        for (int i = data.length - width; i < data.length; i++) {
            smooth[i] = data[i];
        }

        for (int i = width; i < data.length - width; i++) {
            smooth[i] = 0;
            for (int j = 0; j < weights.length; j++) {
                smooth[i] += weights[j] * data[i - width + j];
            }
        }
        return smooth;
    }
}
