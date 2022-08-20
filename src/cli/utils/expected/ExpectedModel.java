package cli.utils.expected;

public abstract class ExpectedModel {

    public abstract double getExpectedFromUncompressedBin(int dist);

    public static double logp1(double x) {
        return Math.log(1 + x);
    }

    public int logp1i(int x) {
        return (int) Math.floor(Math.log(1 + x));
    }

    public double[] expm1(double[] input) {
        double[] vec = new double[input.length];
        for (int k = 0; k < vec.length; k++) {
            vec[k] = Math.expm1(input[k]);
        }
        return vec;
    }
}
