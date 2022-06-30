package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.BitSet;
import java.util.Iterator;


public class Sift {

    private final NormalizationType none = NormalizationHandler.NONE;
    private final NormalizationType scale = NormalizationHandler.SCALE;
    private final int window = 5;

    public Sift(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5);
        }

        // sift <file1.hic> <outfolder>
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1],
                false, false, false);
        File outFolder = UNIXTools.makeDir(new File(args[2]));
        Feature2DList refinedLoops = siftThroughCalls(ds);
        refinedLoops.exportFeatureList(new File(outFolder, "sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private Feature2DList siftThroughCalls(Dataset ds) {
        ChromosomeHandler handler = ds.getChromosomeHandler();

        return null;
    }

    public static void main(String[] args) {
        String file = "/Users/muhammad/Desktop/all_nd.hic";
        Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, false, false);
        //Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr8");
        int maxDist = 10000000;

        for (Chromosome chrom : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            if (chrom.getIndex() < 2 || chrom.getIndex() == 10 || chrom.getIndex() > 20) {
                Matrix matrix = ds.getMatrix(chrom, chrom);
                for (int res : new int[]{200, 100}) {
                    MatrixZoomData zd = matrix.getZoomData(new HiCZoom(res));
                    for (int window : new int[]{5000, 10000}) {
                        int maxBin = maxDist / res;
                        getStats(zd, maxBin, chrom, res, window);
                    }
                }
                matrix.clearCache();
            }
        }

        /*
        int res1 = 100;

        double[] expected = ExpectedUtils.calculateRawExpected(zd, maxBin, true, 0);
        getAllEnrichedPoints(zd, expected);

         */

    }

    private static void getStats(MatrixZoomData zd, int maxBin, Chromosome chrom, int res, int windowSize) {


        //int compression = windowSize / res;
        int maxCompressedBin = logp1i(maxBin) + 1;
        System.out.println("Processing " + chrom.getName() + " " + res + "  max " + maxCompressedBin);

        WelfordStats stats = new WelfordStats(maxCompressedBin);

        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < maxCompressedBin) {
                    stats.addValue(dist, logp1(cr.getCounts()));
                }
            }
        }

        String stem = chrom.getName() + "_" + res + "_w" + windowSize;
        double[] std = stats.getStdDev();
        MatrixTools.saveMatrixTextNumpy(stem + "_std.npy", std);
        QuickMedian.doRollingMedian(std, windowSize / res);
        MatrixTools.saveMatrixTextNumpy(stem + "_smooth_std.npy", std);

        MatrixTools.saveMatrixTextNumpy(stem + "_mu.npy", stats.getMean());
        QuickMedian.doRollingMedian(stats.getMean(), windowSize / res);
        MatrixTools.saveMatrixTextNumpy(stem + "_smooth_mu.npy", stats.getMean());

        int[] positions = new int[maxCompressedBin];
        for (int z = 0; z < maxCompressedBin; z++) {
            positions[z] = (int) Math.expm1(z);
        }

        MatrixTools.saveMatrixTextNumpy(stem + "x.npy", positions);

    }

    private static int logp1i(int x) {
        return (int) Math.floor(Math.log(1 + x));
    }

    private static double logp1(float x) {
        return Math.log(1 + x);
    }


    public static void test1(String[] args) {
        String file = "/Users/muhammad/Desktop/all_nd.hic";
        Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, false, false);
        //Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr8");
        int maxDist = 10000000;

        for (Chromosome chrom : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chrom, chrom);
            for (int res : new int[]{5000, 2000, 1000, 500, 200, 100, 50, 20, 10}) {
                int maxBin = maxDist / res;

                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(res));

                float[] percentages = getPercentNonZeros(zd, maxBin, chrom, res);
                int[] positions = new int[maxBin];
                for (int z = 1; z < maxBin; z++) {
                    positions[z] = positions[z - 1] + res;
                }

                System.out.println("Processing " + chrom.getName() + " " + res);

                String stem = chrom.getName() + "_" + res + "_";
                MatrixTools.saveMatrixTextNumpy(stem + "perc.npy", percentages);
                MatrixTools.saveMatrixTextNumpy(stem + "x.npy", positions);
            }
            matrix.clearCache();
        }

        /*
        int res1 = 100;

        double[] expected = ExpectedUtils.calculateRawExpected(zd, maxBin, true, 0);
        getAllEnrichedPoints(zd, expected);

         */

    }

    private static float[] getPercentNonZeros(MatrixZoomData zd, int maxBin, Chromosome chrom, int res) {

        int maxL = (int) ((chrom.getLength() / res) + 1);
        BitSet rowSums = new BitSet(maxL);

        int[] perc = new int[maxBin];

        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {

                rowSums.set(cr.getBinX());
                rowSums.set(cr.getBinY());

                int dist = ExpectedUtils.getDist(cr);
                if (dist < maxBin) {
                    perc[dist]++;
                }
            }
        }

        int n = rowSums.cardinality();
        rowSums = null;
        float[] percents = new float[maxBin];
        for (int i = 0; i < maxBin; i++) {
            percents[i] = (float) (perc[i]) / (n - i);
        }

        return percents;
    }

    private static void getAllEnrichedPoints(MatrixZoomData zd, double[] expected) {
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {


            }
        }

    }
}
