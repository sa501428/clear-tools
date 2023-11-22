package cli.utils.general.function;

public class ScaledGeoMeanFunction extends NormalizationFunction {
    @Override
    public float[] normalize(float[] original, float[] totals) {
        float[] data = new float[original.length];
        System.arraycopy(original, 0, data, 0, original.length);
        for (int z = 0; z < data.length; z++) {
            if (totals[z] > 0) {
                data[z] = (float) (Math.pow(data[z], 1.0 / totals[z]) * totals[z]);
            }
        }
        return data;
    }

    private double scaleBy(float n, boolean useSqrtScaling) {
        if (useSqrtScaling) {
            return Math.sqrt(n);
        } else {
            return n;
        }
    }
}
