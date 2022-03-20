package cli;

import cli.apa.APA;
import cli.apa.Anchor;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class FLAGS {
    public static void run(String[] args, int resolution, int cutoff, String normString) {

        if(args.length != 4){
            Main.printGeneralUsageAndExit();
        }

        String hicfile = args[1];
        String bedfile = args[2];
        String outfolder = args[3];

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicfile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedfile, cutoff);
        System.out.println("Number of anchors: " + anchors.size());
        APA apa = new APA(ds, outfolder, norm, anchors, resolution);
        apa.run();
    }
}
