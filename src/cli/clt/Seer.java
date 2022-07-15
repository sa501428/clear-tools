package cli.clt;

import cli.utils.seer.CumulativeDistributionFunction;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    // run file with chromosome file, create main class, output as numpy (desktop)

    public static void generateNewReads(String filename, int lowResolution, int highResolution,
                                        String outFolderPath, String possibleNorm) {

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, true);

        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);

        // iterate over a chromosome for now (chromosome 10)
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(highResolution));
            if (zd == null) continue;


            // int[] rowSummation = SeerUtils.getRowSumsForZD(chromosome, highResolution, zd.getDirectIterator());

            MatrixZoomData zdLow = matrix.getZoomData(new HiCZoom(lowResolution));
            if (zdLow == null) continue;

            // create your pdf, cdf (delete the pdf)
            CumulativeDistributionFunction cdf = new CumulativeDistributionFunction(zdLow.getNormalizedIterator(norm),
                    10000000, lowResolution);

            // todo @Allen
            // generate points at random
            // chromosome.getName()
            // <chr1> <pos1> <chr2> <pos2>
            // position =


            //chromToRowSumsMap.put(chromosome, rowSummation);
        }
    }

    public static void run(String[] args, CommandLineParser parser) {
        // check length of arguments equal to 3

        //Map<Chromosome, int[]> chromToRowSumsMap = new HashMap<>();

        int highResolution = parser.getResolutionOption(50);
        int lowResolution = parser.getLowResolutionOption(5000);

        String possibleNorm = parser.getNormalizationStringOption();

        UNIXTools.makeDir(args[2]);

        generateNewReads(args[1], lowResolution, highResolution, args[2], possibleNorm);
        //SeerUtils.exportRowSumsToBedgraph(chromToRowSumsMap, args[2], highResolution);
    }
}