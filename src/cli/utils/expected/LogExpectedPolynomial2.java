package cli.utils.expected;

import cli.utils.sift.ExtremePixels;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class LogExpectedPolynomial2 extends ExpectedModel {

    private static final int minValsPerBin = 10;
    public static int degree = 3;
    private final PolynomialFunction function;
    private final double nearDiagonalSignal;
    private final float maxX;

    public LogExpectedPolynomial2(MatrixZoomData zd, NormalizationType norm,
                                  int maxBin) {
        double[] finalVal = new double[1];
        function = fitDataToFunction(zd, norm, maxBin, finalVal);

        nearDiagonalSignal = Math.expm1(function.value(0.5));
        maxX = findAsymptotePoint(logp1i(maxBin) - 1, finalVal[0]);
    }

    private PolynomialFunction fitDataToFunction(MatrixZoomData zd, NormalizationType norm, int maxBin, double[] finalVal) {

        double[] averageForBin = getAverageInEachBin(zd, norm, maxBin, finalVal);

        List<WeightedObservedPoint> points = new ArrayList<>(averageForBin.length);
        for (int dist0 = 0; dist0 < averageForBin.length; dist0++) {
            if (averageForBin[dist0] > 0) {
                double dist = logp1(dist0);
                double weight = Math.sqrt(1 / (1 + dist));
                if (dist0 < 10) {
                    weight = 10 / (1 + dist);
                }

                points.add(new WeightedObservedPoint(weight, dist, averageForBin[dist0]));
            }
        }
        averageForBin = null;

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);
        return new PolynomialFunction(fitter.fit(points));
    }

    private double[] getAverageInEachBin(MatrixZoomData zd, NormalizationType norm, int maxBin, double[] finalVal) {

        double[] initExpected = new double[maxBin];
        long[] countsPerBin = new long[maxBin];

        float cutoff = (0.75f) * (maxBin);
        double sum = 0;
        long count = 0;

        Iterator<ContactRecord> records = ExtremePixels.getIterator(zd, norm);
        while (records.hasNext()) {
            ContactRecord record = records.next();
            int dist = ExpectedUtils.getDist(record);
            double val = logp1(record.getCounts());
            if (dist < maxBin) {
                initExpected[dist] += val;
                countsPerBin[dist]++;
            }
            if (dist > cutoff) {
                sum += val;
                count++;
            }
        }

        finalVal[0] = sum / count;

        for (int k = 0; k < maxBin; k++) {
            if (countsPerBin[k] > minValsPerBin) {
                initExpected[k] /= countsPerBin[k];
            } else {
                initExpected[k] = -1;
            }
        }

        return initExpected;
    }

    private float findAsymptotePoint(int maxX0, double target) {
        float dx = 0.01f;
        float breakPoint = maxX0;
        for (float x = 1; x < maxX0; x += dx) {
            if (function.value(x - dx) < function.value(x) || function.value(x) < target) {
                breakPoint = x;
                break;
            }
        }
        return breakPoint;
    }

    @Override
    public double getExpectedFromUncompressedBin(int dist0) {
        double dist = Math.min(logp1(dist0), maxX);
        return Math.expm1(function.value(dist));
    }

    @Override
    public double getNearDiagonalSignal() {
        return nearDiagonalSignal;
    }
}
