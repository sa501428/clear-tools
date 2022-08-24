package cli.utils.sample;

import cli.utils.expected.LinearBinnedSpline;
import cli.utils.expected.LogBinnedExpectedModel;
import cli.utils.expected.LogBinnedSpline;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

public class Test2 {

    public static void main(String[] args) {


        Dataset ds = HiCFileTools.extractDatasetForCLT(
                "/Users/muhammad/Desktop/hicfiles/chr10_subsample0.25.hic",
                //"",
                false, false, false);
        Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr10");

        NormalizationType[] norms = new NormalizationType[]{NormalizationHandler.SCALE,
                // NormalizationHandler.VC
        };

        for (int res : new int[]{50000, 10000, 5000, 2000, 1000, 500, 200, 100, 50}) {
            MatrixZoomData zd = ds.getMatrix(chrom, chrom).getZoomData(new HiCZoom(res));

            for (NormalizationType norm : norms) {

                int maxBin = (int) (chrom.getLength() / res);

                LogBinnedExpectedModel model1 = new LogBinnedExpectedModel(zd, norm, maxBin, 0);
                LogBinnedSpline spline = model1.getSpline();
                //LogExpectedPolynomial polynomial = new LogExpectedPolynomial(zd, norm, maxBin);
                LinearBinnedSpline polynomial2 = new LinearBinnedSpline(zd, norm, maxBin, 2);
                LinearBinnedSpline polynomial10 = new LinearBinnedSpline(zd, norm, maxBin, 10);

                System.out.println("Got all points");
                /*
                int n = points.size();
                double[][] data = new double[4][n];
                int counter = 0;
                for (WeightedObservedPoint point : points) {
                    data[0][counter] = point.getX();
                    data[1][counter] = point.getY();
                    counter++;
                }
                */

                double[][] data = new double[6][maxBin];
                double[] evec = ds.getExpectedValues(new HiCZoom(res), norm, false).getExpectedValuesWithNormalization(chrom.getIndex()).getValues().get(0);
                for (int x = 0; x < maxBin; x++) {
                    data[0][x] = x;
                    data[1][x] = model1.getExpectedFromUncompressedBin(x);
                    data[2][x] = spline.getExpectedFromUncompressedBin(x);
                    data[3][x] = polynomial2.getExpectedFromUncompressedBin(x);
                    data[4][x] = polynomial10.getExpectedFromUncompressedBin(x);
                    data[5][x] = evec[x];
                }
                MatrixTools.saveMatrixTextNumpy("multi_interp_" + norm.getLabel() + "_" + res + ".npy", data);
            }
        }
    }
}
