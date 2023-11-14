package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.anchors.AnchorFrequencyCounter;
import cli.utils.anchors.Convolution1DTools;
import cli.utils.anchors.LoopToAnchorIntersecter;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.awt.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class IntegrateLoopListsAndUnWrap {
    private static final int numLists = 6;
    private static final int lowResolution = 1000;
    private static final int searchWidth = 250;
    // [-filter]
    public static String usage = "integrate-loops [-r hi-res] <genomeID> <output.bedpe> <file1.bedpe> [file2.bedpe] ... [fileN.bedpe]\n" +
            "\t\tunwrap localizer output to proper hi-res inverted bounds list\n" +
            "\t\tlegacy mode will use the old method of unwrapping the combined anchors\n";
    private static int highResolution = 10;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length < 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        highResolution = parser.getResolutionOption(10);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        String outFile = args[2];

        List<Feature2DList> loopLists = new LinkedList<>();
        for (int z = 3; z < args.length; z++) {
            loopLists.add(Feature2DParser.loadFeatures(args[z], handler,
                    true, null, false));
        }

        Map<String, int[]> chromToHiResAnchorCounts = AnchorFrequencyCounter.calculateHiResAnchorCounts(loopLists,
                handler, highResolution);

        Map<String, float[]> chromToHiResAnchorSmoothCounts = smoothCounts(chromToHiResAnchorCounts);

        Feature2DList mergedDeDuped = LoopToAnchorIntersecter.intersect(loopLists, handler,
                chromToHiResAnchorSmoothCounts, highResolution);
    }

    private static Map<String, float[]> smoothCounts(Map<String, int[]> chromToHiResAnchorCounts) {
        Map<String, float[]> chromToHiResAnchorSmoothCounts = new HashMap<>();
        for (String key : chromToHiResAnchorCounts.keySet()) {
            int[] counts = chromToHiResAnchorCounts.get(key);
            chromToHiResAnchorSmoothCounts.put(key, Convolution1DTools.smooth5(counts));
        }
        return chromToHiResAnchorSmoothCounts;
    }


    private static Feature2D unwrapLocal(Feature2D feature2D) {
        try {
            long start1 = Long.parseLong(feature2D.getAttribute("localX"));
            long start2 = Long.parseLong(feature2D.getAttribute("localY"));
            long end1 = start1 + highResolution;
            long end2 = start2 + highResolution;
            Map<String, String> attrs = new HashMap<>(feature2D.getAttributes());
            return new Feature2D(Feature2D.FeatureType.PEAK, feature2D.getChr1(), start1, end1,
                    feature2D.getChr2(), start2, end2, Color.BLUE, attrs);
        } catch (Exception ignored) {
        }
        return null;
    }

}
