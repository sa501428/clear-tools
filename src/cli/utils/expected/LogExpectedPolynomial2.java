package cli.utils.expected;

import cli.utils.sift.ExtremePixels;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class LogExpectedPolynomial2 extends ExpectedModel {

    private static final int minValsPerBin = 10;
    public static int degree = 3;
    private final PolynomialSplineFunction function;
    private final double nearDiagonalSignal;
    private final float maxX;

    public LogExpectedPolynomial2(MatrixZoomData zd, NormalizationType norm,
                                  int maxBin, int base) {
        double[] finalVal = new double[1];
        function = fitDataToFunction(zd, norm, maxBin, finalVal, base);
        nearDiagonalSignal = Math.expm1(function.value(0.5));
        maxX = logp1i(maxBin);
        //maxX = findAsymptotePoint(logp1i(maxBin) - 1, finalVal[0]);
    }

    private PolynomialSplineFunction fitDataToFunction(MatrixZoomData zd, NormalizationType norm, int maxBin,
                                                       double[] finalVal, int base) {

        double[] averageForBin = getAverageInEachBin(zd, norm, maxBin, finalVal);

        List<double[]> points = new ArrayList<>(averageForBin.length / 2);
        int dx = 1;
        int counter = 0;
        for (int dist0 = 0; dist0 < averageForBin.length; dist0 += dx) {
            double dist = logp1(dist0);
            if (dx > 1) {
                double newVal = QuickMedian.fastMedian(extract(averageForBin, dist0, dx));
                if (newVal > 0) {
                    points.add(new double[]{dist, newVal});
                }
            } else {
                if (averageForBin[dist0] > 0) {
                    points.add(new double[]{dist, averageForBin[dist0]});
                }
            }
            if (++counter % base == 0) {
                dx *= base;
            }
        }
        points.add(new double[]{logp1(maxBin), finalVal[0]});
        points.add(new double[]{logp1(maxBin) + 1, finalVal[0]});
        averageForBin = null;

        double[] x = new double[points.size()];
        double[] y = new double[points.size()];
        for (int i = 0; i < points.size(); i++) {
            x[i] = points.get(i)[0];
            y[i] = points.get(i)[1];
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        return interpolator.interpolate(x, y);
    }

    private List<Double> extract(double[] array, int x0, int dx) {
        List<Double> positives = new ArrayList<>(dx * 2 + 1);
        for (int x = Math.max(0, x0 - dx); x < Math.min(x0 + dx + 1, array.length); x++) {
            if (array[x] > 0) {
                positives.add(array[x]);
            }
        }
        return positives;
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
        double dist = Math.max(0, logp1(dist0));
        dist = Math.min(dist, maxX);
        return Math.expm1(function.value(dist));
    }

    @Override
    public double getNearDiagonalSignal() {
        return nearDiagonalSignal;
    }
}
