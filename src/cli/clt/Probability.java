package cli.clt;

import cli.Main;
import cli.utils.general.HiCValue;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.List;

public class Probability {

    // "probability [--res int] <input.hic> <loops.bedpe> <out_folder>\n" +

    private static final int DEFAULT_MAX_DIST = 20000000;
    private static final int DEFAULT_WINDOW = 1000000;


    public static void run(String[] args, CommandLineParser parser) {

        if (args.length != 4) {
            Main.printGeneralUsageAndExit(4, null);
        }

        String hicFile = args[1];
        String bedpeFile = args[2];
        String outFolder = args[3];

        int resolution = parser.getResolutionOption(5000);
        boolean useLog = parser.getLogOption();

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicFile, false, true, true);
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

            double[] expected = ExpectedUtils.calculateExpected(zd, norm, DEFAULT_MAX_DIST / resolution, useLog);
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
            int dist = ExpectedUtils.getDist(loop, resolution);
            float val = HiCValue.getExactPixel(zd, loop, resolution, norm);
            double proportionLoop = (val - expected[dist])/(expected[1]-expected[dist]);
            double enrichment = val / expected[dist];
            double pseudocount = expected[expected.length-1];
            //double enrichment2 = (val + 1) / (expected[dist] + 1);
            double enrichment3 = (val + pseudocount) / (expected[dist] + pseudocount);

            System.out.println("Loop "+i+" "+loop.simpleString()+" pL "+proportionLoop +
                    " O/E "+enrichment +
                    //" O+1/E+1 "+enrichment2 +
                    " O+ps/E+ps "+enrichment3 +
                    " ps "+pseudocount
            );
        }
    }
}
