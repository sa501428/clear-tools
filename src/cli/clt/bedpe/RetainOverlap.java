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

public class RetainOverlap {
    public static String usage = "retain-exact-overlap <genomeID> <file1.bedpe> <file2.bedpe> " +
            "<output.bedpe>\n" +
            "behavior retains loops (with attributes) in file1 that overlap exactly with" +
            "a loop in file2";

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(58, usage);
        }

        String genomeID = args[1];
        String inputPath1 = args[2];
        String inputPath2 = args[3];
        String outPath = args[4];

        simpleRetain(inputPath1, inputPath2, genomeID, outPath);

        System.out.println("retain complete");
    }

    private static void simpleRetain(String inputPath1, String inputPath2, String genomeID, String outPath) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList1 = Feature2DParser.loadFeatures(inputPath1, handler,
                true, null, false);
        Feature2DList loopList2 = Feature2DParser.loadFeatures(inputPath2, handler,
                false, null, false);

        Feature2DList retained = new Feature2DList();
        loopList1.processLists((key, list) -> {
            Set<String> anchors = getLoopCodes(loopList2.get(key));
            List<Feature2D> toKeep = new LinkedList<>();
            for (Feature2D loop : list) {
                if (anchors.contains(getLoopCode(loop))) {
                    toKeep.add(loop);
                }
            }
            retained.addByKey(key, toKeep);
        });

        retained.exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);
    }

    private static Set<String> getLoopCodes(List<Feature2D> feature2DS) {
        Set<String> codes = new HashSet<>();
        for (Feature2D loop : feature2DS) {
            codes.add(getLoopCode(loop));
        }
        return codes;
    }

    private static String getLoopCode(Feature2D loop) {
        return loop.getChr1() + "_" + loop.getStart1() + "_" + loop.getEnd1() + "_" +
                loop.getChr2() + "_" + loop.getStart2() + "_" + loop.getEnd2();
    }

}
