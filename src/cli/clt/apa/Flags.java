package cli.clt.apa;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.flags.Anchor;
import cli.utils.flags.FlagsAggregation;
import cli.utils.general.BedFileParser;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class Flags {
    public static String usage = "flags [--cutoff int] [--res int] [--norm string] <input.hic> <loops.bedpe> <out_folder>";

    public static void run(String[] args, CommandLineParser parser) {

        if (args.length != 4) {
            Main.printGeneralUsageAndExit(4, usage);
        }

        String hicFile = args[1];
        String bedFile = args[2];
        String outFolder = args[3];

        int resolution = parser.getResolutionOption(5000);
        int cutoff = parser.getCutoffOption(500);
        String normString = parser.getNormalizationStringOption();

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicFile, false, true, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        NormalizationType norm;
        if (normString != null && normString.length() > 0) {
            try {
                norm = ds.getNormalizationHandler().getNormTypeFromString(normString);
            } catch (Exception e) {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "SCALE", "KR", "NONE"});
            }
        } else {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "NONE"});
        }
        System.out.println("Norm being used: " + norm.getLabel());

        GenomeWide1DList<Anchor> anchors = BedFileParser.loadFromBEDFile(handler, bedFile, cutoff, false);
        System.out.println("Number of anchors: " + anchors.size());
        FlagsAggregation apa = new FlagsAggregation(ds, outFolder, norm, anchors, resolution);
        apa.run();
    }
}
