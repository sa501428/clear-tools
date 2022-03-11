package flags;

import flags.apa.APA;
import flags.apa.Anchor;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class AVD {

    public AVD(String hicfile, String bedfile) {


        Dataset ds = HiCFileTools.extractDatasetForCLT(hicfile, false, false);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"KR", "SCALE", "NONE"});

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedfile, 500);
        APA apa = new APA(ds, "apa_result", norm, anchors);
        apa.run();
    }
}
