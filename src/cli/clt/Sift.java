package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.sift.SimpleLocation;
import cli.utils.sift.ZScores;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;


public class Sift {
    private final int MAX_DIST = 10000000;
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

    public static void main(String[] args) {
        String file = "/Users/muhammad/Desktop/all_nd.hic";
        //Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr8");
        int maxDist = 10000000;



        /*
        int res1 = 100;

        double[] expected = ExpectedUtils.calculateRawExpected(zd, maxBin, true, 0);
        getAllEnrichedPoints(zd, expected);

         */

    }

    private static Set<SimpleLocation> getExtremePixels(MatrixZoomData zd, int maxBin) {

        int maxCompressedBin = logp1i(maxBin) + 1;
        ZScores zScores = getZscores(zd, maxCompressedBin);
        Set<SimpleLocation> records = new HashSet<>();

        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < maxCompressedBin) {
                    float zscore = zScores.getZscore(dist, logp1(cr.getCounts()));
                    if (zscore > 3) {
                        records.add(new SimpleLocation(cr));
                    }
                }
            }
        }

        return records;
    }

    private static ZScores getZscores(MatrixZoomData zd, int length) {
        WelfordStats stats = new WelfordStats(length);
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < length) {
                    stats.addValue(dist, logp1(cr.getCounts()));
                }
            }
        }
        return stats.getZscores();
    }

    private Feature2DList siftThroughCalls(Dataset ds) {
        ChromosomeHandler handler = ds.getChromosomeHandler();

        for (Chromosome chrom : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            if (chrom.getIndex() < 2 || chrom.getIndex() == 10 || chrom.getIndex() > 20) {
                Matrix matrix = ds.getMatrix(chrom, chrom);

                MatrixZoomData zd1 = matrix.getZoomData(new HiCZoom(200));
                Set<SimpleLocation> initialPoints = getExtremePixels(zd1, MAX_DIST / 200);

                MatrixZoomData zd2 = matrix.getZoomData(new HiCZoom(100));
                Set<SimpleLocation> finalPoints = getOverlappingExtremePixels(zd2,
                        MAX_DIST / 100, initialPoints, 200 / 100);

                initialPoints.clear();

                matrix.clearCache();
            }
        }

        return null;
    }

    private Set<SimpleLocation> getOverlappingExtremePixels(MatrixZoomData zd, int maxBin, Set<SimpleLocation> regions, int scalar) {

        int maxCompressedBin = logp1i(maxBin) + 1;
        ZScores zScores = getZscores(zd, maxCompressedBin);
        Set<SimpleLocation> records = new HashSet<>();

        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < maxCompressedBin) {
                    float zscore = zScores.getZscore(dist, logp1(cr.getCounts()));
                    if (zscore > 3) {
                        if (inRegions(cr, regions, scalar)) {
                            records.add(new SimpleLocation(cr));
                        }
                    }
                }
            }
        }

        return records;
    }

    private boolean inRegions(ContactRecord cr, Set<SimpleLocation> regions, int scalar) {
        SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
        return regions.contains(region);
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
