package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Clique {

    // [-rescue]
    public static String usage = "clique[-clean][-rescue] [-r resolution] <genomeID> <input.bedpe> <output.bedpe>\n" +
            "\t\tdefault behavior finds the cliques using midpoints of the anchors at the resolution specified\n" +
            "\t\trescue will predict loops that were potentially missed\n" +
            "\t\tclean avoids saving old attributes";
    private static int resolution = 500;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }
        resolution = parser.getResolutionOption(500);
        String genomeID = args[1];
        String inFile = args[2];
        String outStem = args[3];
        getCliques(inFile, genomeID, outStem,
                command.contains("clean"), command.contains("rescue"));
        System.out.println("anchor expansion complete");
    }

    private static void getCliques(String inputBedpe, String genomeID, String outStem,
                                   boolean noAttributes, boolean rescueLoops) {
        if (rescueLoops) {
            System.err.println("Not yet implemented");
            System.exit(9);
        }

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpe, handler, !noAttributes,
                null, false);

        loopList.filterLists((s, list) -> {
            return recoverLoops(list);
        });

        loopList.exportFeatureList(new File(outStem + ".rescue.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> recoverLoops(List<Feature2D> list) {
        List<Feature2D> newLoops = new ArrayList<>(2 * list.size());

        Map<Integer, Integer> upStreamAnchors = new HashMap<>();
        Map<Integer, Integer> downStreamAnchors = new HashMap<>();


        // matrix multiplication
        // N^3

        /*
        for(int i = 0; i < n; i++){
            for(int j = i+1; j < n; j++){

            }
        }
        */

        return newLoops;
    }


    private static Feature2D getUpdatedLoop(Feature2D feature, int newAnchorSize, boolean setExactSize) {
        long mid1 = feature.getMidPt1();
        long mid2 = feature.getMidPt2();
        long start1 = mid1 - (newAnchorSize / 2);
        long start2 = mid2 - (newAnchorSize / 2);
        long end1 = start1 + newAnchorSize;
        long end2 = start2 + newAnchorSize;

        if (!setExactSize) {
            if (feature.getWidth1() > newAnchorSize) {
                start1 = feature.getStart1();
                end1 = feature.getEnd1();
            }
            if (feature.getWidth2() > newAnchorSize) {
                start2 = feature.getStart2();
                end2 = feature.getEnd2();
            }
        }
        return new Feature2D(feature.getFeatureType(),
                feature.getChr1(), start1, end1,
                feature.getChr2(), start2, end2,
                feature.getColor(), feature.getAttributes());
    }
}
