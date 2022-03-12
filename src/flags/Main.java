package flags;

import flags.apa.APA;
import flags.apa.Anchor;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class Main {

    public static void main(String[] args) {
        if (args.length < 6) {
            System.err.println("Usage: flags <hicfile> <bedfile> <score_cutoff> <resolution> <norm> <outfolder>");
            System.exit(8);
        }

        String hicfile = args[0];
        String bedfile = args[1];
        int cutoff = Integer.parseInt(args[2]);
        int resolution = Integer.parseInt(args[3]);
        String normString = args[4];
        String outfolder = args[5];

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
