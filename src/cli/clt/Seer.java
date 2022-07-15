package cli.clt;

import cli.utils.seer.SeerUtils;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    // run file with chromosome file, create main class, output as numpy (desktop)

    public static void calculateRowSums(String filename, Map<Chromosome, int[]> chromToRowSumsMap, int lowResolution, int highResolution) {

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, true);

        // iterate over a chromosome for now (chromosome 10)
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(highResolution));
            if (zd == null) continue;

            int[] rowSummation = SeerUtils.getRowSumsForZD(chromosome, highResolution, zd.getDirectIterator());

            chromToRowSumsMap.put(chromosome, rowSummation);
        }
    }

    public static void run(String[] args, CommandLineParser parser) {
        // check length of arguments equal to 3

        Map<Chromosome, int[]> chromToRowSumsMap = new HashMap<>();

        int highResolution = parser.getResolutionOption(50);
        int lowResolution = parser.getLowResolutionOption(5000);

        calculateRowSums(args[1], chromToRowSumsMap, lowResolution, highResolution);
        UNIXTools.makeDir(args[2]);
        try {
            SeerUtils.exportRowSumsToBedgraph(chromToRowSumsMap, args[2], highResolution);
        } catch (IOException e) {
            e.printStackTrace();
            System.err.println("Unable to export sums to bedgraph");
        }
    }
}