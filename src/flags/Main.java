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
        String hicfile = args[0];
        String bedfile = args[1];
        int cutoff = Integer.parseInt(args[2]);
        int resolution = Integer.parseInt(args[3]);
        String normString = args[4];

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicfile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "KR", "SCALE", "NONE"});

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedfile, cutoff);
        System.out.println("Number of anchors: " + anchors.size());
        APA apa = new APA(ds, "apa_result", norm, anchors, resolution);
        apa.run();
    }
}
