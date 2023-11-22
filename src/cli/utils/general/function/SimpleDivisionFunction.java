package cli.utils.general.function;

public class SimpleDivisionFunction extends NormalizationFunction {
    @Override
    public float[] normalize(float[] data, float[] totals) {
        float[] normalized = new float[data.length];
        for (int i = 0; i < data.length; i++) {
            normalized[i] = data[i] / totals[i];
        }
        return normalized;
    }
}
