package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class SubtractByAnchorOverlap {

    public static String usage = "subtract-by-anchor-overlap [-r resolution] <genomeID> <A.bedpe> <B.bedpe> <output.stem>\n" +
            "\t\tdefault behavior removes any loop from A that overlaps anchors " +
            "with a loop in B.\n";
    private static int resolution = 5000;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(15, usage);
        }
        String genomeID = args[1];
        String inFileA = args[2];
        String inFileB = args[3];
        String outStem = args[4];
        resolution = parser.getResolutionOption(resolution);
        filter(inFileA, inFileB, genomeID, outStem);
        System.out.println("filtering complete");
    }

    private static void filter(String inputBedpeA, String inputBedpeB, String genomeID, String outStem) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopListA = Feature2DParser.loadFeatures(inputBedpeA, handler, true,
                null, false);
        Feature2DList loopListB = Feature2DParser.loadFeatures(inputBedpeB, handler, true,
                null, false);

        Feature2DList aMinusB = new Feature2DList();
        Feature2DList whatWasRemoved = new Feature2DList();

        loopListA.processLists((key, list) -> {
            Set<String> anchors = getAnchorSet(loopListB.get(key));

            List<Feature2D> toKeep = new LinkedList<>();
            List<Feature2D> toRemove = new LinkedList<>();
            for (Feature2D loop : list) {
                if (anchors.contains(getAnchorKey(loop.getChr1(), loop.getStart1(), loop.getEnd1())) &&
                        anchors.contains(getAnchorKey(loop.getChr2(), loop.getStart2(), loop.getEnd2()))) {
                    toRemove.add(loop);
                } else {
                    toKeep.add(loop);
                }
            }

            aMinusB.addByKey(key, toKeep);
            whatWasRemoved.addByKey(key, toRemove);
        });

        aMinusB.exportFeatureList(new File(outStem + ".bedpe"), false, Feature2DList.ListFormat.NA);
        whatWasRemoved.exportFeatureList(new File(outStem + ".removed.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static Set<String> getAnchorSet(List<Feature2D> loops) {
        Set<String> anchors = new HashSet<>();
        for (Feature2D loop : loops) {
            anchors.add(getAnchorKey(loop.getChr1(), loop.getStart1(), loop.getEnd1()));
            anchors.add(getAnchorKey(loop.getChr2(), loop.getStart2(), loop.getEnd2()));
        }
        return anchors;
    }

    private static String getAnchorKey(String chr, long start, long end) {
        int mid = (int) (((start + end) / 2) / resolution);
        return chr + "_" + mid;
    }
}
