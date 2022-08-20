package cli.utils.expected;

import cli.utils.sift.ExtremePixels;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

public class LogExpectedPolynomial extends ExpectedModel {

    private final PolynomialFunction function;
    private final double nearDiagonalSignal;
    private final float maxX;

    public LogExpectedPolynomial(MatrixZoomData zd, NormalizationType norm,
                                 int maxBin, int degree, boolean useAsymptote) {

        Random generator = new Random(0);
        Iterator<ContactRecord> records = ExtremePixels.getIterator(zd, norm);
        List<WeightedObservedPoint> points = new LinkedList<>();

        while (records.hasNext()) {
            ContactRecord record = records.next();
            if (generator.nextDouble() < 0.01) {
                double dist = logp1(ExpectedUtils.getDist(record));
                double val = logp1(record.getCounts());
                points.add(new WeightedObservedPoint(Math.sqrt(1.0 / (1.0 + dist)), dist, val)); //
            }
        }


        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);
        function = new PolynomialFunction(fitter.fit(points));
        nearDiagonalSignal = Math.expm1(function.value(0.5));

        if (useAsymptote) {
            maxX = findAsymptotePoint(logp1i(maxBin) - 1);
        } else {
            maxX = logp1i(maxBin) - 1;
        }
    }

    private float findAsymptotePoint(int maxX0) {
        float dx = 0.01f;
        float breakPoint = maxX0;
        for (float x = 1; x < maxX0; x += dx) {
            if (function.value(x - dx) < function.value(x)) {
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
