package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class AnchorFix {

    private static final int widthToConnect = 2;
    public static String usage = "anchor-fix[-clean] [-r resolution] <genomeID> <input.bedpe> <output.stem>\n" +
            "\t\tdefault behavior will fix the hi-res shared anchors for loops\n" +
            "\t\tclean avoids saving old attributes";
    private static int resolution = 100;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(57, usage);
        }
        resolution = parser.getResolutionOption(resolution);
        String genomeID = args[1];
        String inFile = args[2];
        String outStem = args[3];
        fixAnchors(inFile, genomeID, outStem, command.contains("clean"));
    }


    private static void fixAnchors(String inputBedpe, String genomeID, String outStem, boolean noAttributes) {
        Feature2DList output = new Feature2DList();
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpe, handler, !noAttributes,
                null, false);
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
            List<Feature2D> loops = loopList.get(chrom.getIndex(), chrom.getIndex());
            if (loops.size() > 0) {
                if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
                List<Feature2D> newLoops = recoverLoops(loops);
                output.addByKey(Feature2DList.getKey(chrom, chrom), newLoops);
            }
        }

        output.exportFeatureList(new File(outStem + ".fixed.anchors.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> recoverLoops(List<Feature2D> initialLoops) {
        return new ArrayList<>();
    }
}
