package cli.utils.general.function;

import cli.utils.general.VectorCleaner;

public class OELogMedianFunction extends NormalizationFunction {
    @Override
    public float[] normalize(float[] data, float[] totals) {
        return VectorCleaner.logMedianScale(data);
    }
}
