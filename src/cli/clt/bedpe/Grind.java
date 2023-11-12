package cli.clt.bedpe;

import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.io.File;

public class Grind {

    public static String usage = "grind [-k NORM] [-r resolution] [--window half-width] <hic file> <bedpe> <directory>\n" +
            "\t\tsplit up the bedpe into multiple lists; number < 2 splits by chromosome";
    private final boolean useObservedOverExpected = false;
    private final Dataset ds;
    private final File outputDirectory;
    private final HiCZoom zoom;
    private final Feature2DList loopList;
    private final ChromosomeHandler handler;
    private int matrixHalfWidth = 10;
    private Integer resolution = 1000;
    private NormalizationType norm;

    public Grind(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(5000);
        boolean useBI = resolution >= 50;

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, useBI);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        String possibleNorm = parser.getNormalizationStringOption();

        try {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
        }

        System.out.println("Using normalization: " + norm.getLabel());


        matrixHalfWidth = parser.getWindowSizeOption(10);

        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();


        loopList = Feature2DParser.loadFeatures(args[2], handler, false, null, false);
        if (loopList.getNumTotalFeatures() < 1) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }
    }

    public void run() {
        LoopDumper.dump(ds, loopList, outputDirectory, handler, norm,
                useObservedOverExpected, resolution, matrixHalfWidth);
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }
}
