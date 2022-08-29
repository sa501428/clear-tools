package cli.utils.sample;

import javastraw.expected.ExpectedUtils;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.util.BitSet;
import java.util.Iterator;

public class Test1 {

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
}
