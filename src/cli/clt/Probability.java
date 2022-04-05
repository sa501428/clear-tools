package cli.clt;

import cli.HiCValue;
import cli.Main;
import cli.apa.APAUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.Iterator;
import java.util.List;

public class Probability {
    private static final int DEFAULT_MAX_DIST = 20000000;
    private static final int DEFAULT_WINDOW = 1000000;


    public static void run(String[] args, int resolution, boolean useLog) {

        if(args.length != 4){
            Main.printGeneralUsageAndExit(4);
        }

        String hicFile = args[1];
        String bedpeFile = args[2];
        String outFolder = args[3];

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicFile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(bedpeFile, handler, false, null, false);
        UNIXTools.makeDir(outFolder);

        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        calculateProbabilities(ds, handler, loopList, outFolder, norm, resolution, useLog);
    }

    private static void calculateProbabilities(Dataset ds, ChromosomeHandler handler,
                                               Feature2DList loopList, String outFolder,
                                               NormalizationType norm, int resolution, boolean useLog) {

        for(Chromosome chromosome : handler.getChromosomeArrayWithoutAllByAll()) {
            if (chromosome.getIndex() != 10) continue; // temporary
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) continue;

            double[] expected = calculateExpected(zd, norm, DEFAULT_MAX_DIST/resolution, useLog);
            MatrixTools.saveMatrixTextNumpy((new File(outFolder, "expected.npy")).getAbsolutePath(),
                    expected);
            QuickMedian.doRollingMedian(expected, DEFAULT_WINDOW/resolution);
            MatrixTools.saveMatrixTextNumpy((new File(outFolder, "smooth_expected.npy")).getAbsolutePath(),
                    expected);
            getPercentages(expected, loopList.get(chromosome.getIndex(), chromosome.getIndex()),
                    zd, resolution, norm);
        }

    }

    private static void getPercentages(double[] expected, List<Feature2D> loops, MatrixZoomData zd,
                                       int resolution, NormalizationType norm) {
        for(int i = 0; i < 5; i++) {
            Feature2D loop = loops.get(i);
            int dist = getDist(loop, resolution);
            float val = HiCValue.getExactPixel(zd, loop, resolution, norm);
            double proportionLoop = (val - expected[dist])/(expected[1]-expected[dist]);
            double enrichment = val / expected[dist];
            double pseudocount = expected[expected.length-1];
            //double enrichment2 = (val + 1) / (expected[dist] + 1);
            double enrichment3 = (val + pseudocount) / (expected[dist] + pseudocount);
            /*
            System.out.println("Loop "+i+" "+loop.simpleString()+" pL "+proportionLoop +
                    " O/E "+enrichment +
                    //" O+1/E+1 "+enrichment2 +
                    " O+ps/E+ps "+enrichment3 +
                    " ps "+pseudocount
            );
            */
            System.out.println("Loop "+i+" "+loop.simpleString()+" pL "+proportionLoop +
                    " O/E "+enrichment +
                    //" O+1/E+1 "+enrichment2 +
                    " O+ps/E+ps "+enrichment3 +
                    " ps "+pseudocount
            );
        }
    }

    private static double[] calculateExpected(MatrixZoomData zd, NormalizationType norm, int maxBinDist,
                                              boolean useLog) {
        double[] expected = new double[maxBinDist];
        long[] counts = new long[maxBinDist];

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()){
            ContactRecord record = iterator.next();
            int dist = getDist(record);
            if(dist < maxBinDist){
                if(useLog){
                    expected[dist] += Math.log(1 + record.getCounts());
                } else {
                    expected[dist] += record.getCounts();
                }
                counts[dist]++;
            }
        }

        for(int z = 0; z < maxBinDist; z++){
            if(counts[z] > 0){
                expected[z] /= counts[z];
            }
        }

        if(useLog){
            for(int z = 0; z < maxBinDist; z++){
                expected[z] = Math.expm1(expected[z]);
            }
        }

        return expected;
    }

    private static int getDist(Feature2D loop, int resolution) {
        int binXStart = (int) (loop.getMidPt1() / resolution);
        int binYStart = (int) (loop.getMidPt2() / resolution);
        return Math.abs(binXStart-binYStart);
    }

    private static int getDist(ContactRecord record) {
        return Math.abs(record.getBinX()- record.getBinY());
    }
}
