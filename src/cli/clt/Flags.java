package cli.clt;

import cli.utils.GenericLocusParser;
import cli.Main;
import cli.apa.APA;
import cli.apa.Anchor;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class Flags {
    public static void run(String[] args, int resolution, int cutoff, String normString) {

        if(args.length != 4){
            Main.printGeneralUsageAndExit(4);
        }

        String hicFile = args[1];
        String bedFile = args[2];
        String outFolder = args[3];

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicFile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedFile, cutoff);
        System.out.println("Number of anchors: " + anchors.size());
        APA apa = new APA(ds, outFolder, norm, anchors, resolution);
        apa.run();
    }
}
