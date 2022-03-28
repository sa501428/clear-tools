package cli.clt;

import cli.Main;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;

public class Probability {
    public static void run(String[] args, int resolutionOption) {

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

        

    }
}
