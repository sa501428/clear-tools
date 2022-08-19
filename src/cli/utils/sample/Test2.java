package cli.utils.sample;

import cli.utils.expected.ExpectedUtils;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

public class Test2 {

    public static void main(String[] args) {

        int res = 5000;

        Dataset ds = HiCFileTools.extractDatasetForCLT(
                "/Users/muhammad/Desktop/hicfiles/chr10_subsample0.25.hic",
                false, false, false);
        Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr10");
        MatrixZoomData zd = ds.getMatrix(chrom, chrom).getZoomData(new HiCZoom(res));
        Iterator<ContactRecord> records = zd.getNormalizedIterator(NormalizationHandler.SCALE);


        List<WeightedObservedPoint> points = new LinkedList<>();
        Random generator = new Random(0);

        while (records.hasNext()) {
            ContactRecord record = records.next();
            //if(record.getCounts() > 1) {
            int dist0 = ExpectedUtils.getDist(record);
            //if (dist0 < 10000000/res) {
            if (generator.nextDouble() < 0.001) {
                double dist = (Math.log(1 + dist0));
                double val = Math.log(1 + record.getCounts());
                points.add(new WeightedObservedPoint(Math.sqrt(1.0 / (1.0 + dist)), dist, val)); //
            }
            //}
            //}
        }

        System.out.println("Got all points");

        int n = points.size();
        double[][] data = new double[4][n];

        int counter = 0;
        for (WeightedObservedPoint point : points) {
            data[0][counter] = point.getX();
            data[1][counter] = point.getY();
            counter++;
        }

        for (int k = 0; k < n; k++) {
            double x = ((double) k * Math.log(100000000 / res)) / ((double) n - 1);
            data[2][k] = x;
        }

        data[3] = fitPolynomial(points, data[2], 5);

        MatrixTools.saveMatrixTextNumpy("temp_interp.npy", data);
    }

    public static double[] fitPolynomial(double[] x0, double[] y0,
                                         double[] xF, int degree) {

        List<WeightedObservedPoint> points = new LinkedList<>();

        for (int i = 0; i < x0.length; i++) {
            points.add(new WeightedObservedPoint(1.0, x0[i], y0[i]));
        }

        return fitPolynomial(points, xF, degree);
    }

    public static double[] fitPolynomial(List<WeightedObservedPoint> points,
                                         double[] xF, int degree) {

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);

        double[] coeff = fitter.fit(points);

        PolynomialFunction function = new PolynomialFunction(coeff);

        double[] fitted = new double[xF.length];

        for (int i = 0; i < xF.length; i++) {
            fitted[i] = function.value(xF[i]);
        }

        return fitted;
    }

    public void main2(String[] args) {
        Random generator = new Random();
        int n = 100000;
        double[][] data = new double[4][n];
        data[0][1] = 5;
        data[1][1] = 25;

        for (int k = 2; k < n; k++) {
            data[0][k] = generator.nextDouble() * 5;
            data[1][k] = (data[0][k] * data[0][k] + (data[0][k] * 2) * generator.nextGaussian() / 2);
        }

        for (int k = 0; k < n; k++) {
            double x = ((double) k * 5.0) / ((double) n - 1);
            data[2][k] = x;
        }

        data[3] = fitPolynomial(data[0], data[1], data[2], 4);

        MatrixTools.saveMatrixTextNumpy("temp_interp.npy", data);
    }
}
