package cli.utils.expected;

import javastraw.reader.block.ContactRecord;

public abstract class ExpectedModel {

    public abstract double getExpectedFromUncompressedBin(int dist);

    public static float getP(double obs, double expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    public abstract double getNearDiagonalSignal();

    public float getPercentContact(ContactRecord cr) {
        return getPercentContact(ExpectedUtils.getDist(cr), cr.getCounts());
    }

    public static double logp1(double x) {
        return Math.log(1 + x);
    }

    public float getPercentContact(int dist, float counts) {
        double baseline = getExpectedFromUncompressedBin(dist);
        return getP(counts, baseline, getNearDiagonalSignal());
    }

    public double[] expm1(double[] input) {
        double[] vec = new double[input.length];
        for (int k = 0; k < vec.length; k++) {
            vec[k] = Math.expm1(input[k]);
        }
        return vec;
    }

    public int logp1i(int x) {
        return (int) Math.log(1 + x);
    }

    public boolean isReasonablePercentContact(ContactRecord cr) {
        double percentContact = getPercentContact(cr);
        return percentContact > 0.01;// && percentContact < 0.4;
    }
}
