package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.anchors.AnchorFrequencyCounter;
import cli.utils.anchors.Convolution1DTools;
import cli.utils.anchors.LoopToAnchorIntersecter;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class IntegrateLoopListsAndUnWrap {
    private static final int numLists = 6;
    private static final int lowResolution = 1000;
    private static final int searchWidth = 250;
    // [-filter]
    public static String usage = "integrate-loops [-r hi-res] [--window window] <genomeID> <output.bedpe> <file1.bedpe> [file2.bedpe] ... [fileN.bedpe]\n" +
            "\t\tunwrap localizer output to proper hi-res inverted bounds list\n" +
            "\t\tlegacy mode will use the old method of unwrapping the combined anchors\n";
    private static int highResolution = 10;

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length < 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        highResolution = parser.getResolutionOption(10);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        String outPath = args[2];

        List<Feature2DList> loopLists = new LinkedList<>();
        for (int z = 3; z < args.length; z++) {
            loopLists.add(Feature2DParser.loadFeatures(args[z], handler,
                    true, null, false));
        }

        int window = parser.getWindowSizeOption(10);

        Map<String, int[]> chromToHiResAnchorCounts = AnchorFrequencyCounter.calculateHiResAnchorCounts(loopLists,
                handler, highResolution);

        Map<String, float[]> chromToHiResAnchorSmoothCounts = smoothCounts(chromToHiResAnchorCounts);

        Feature2DList mergedDeDuped = LoopToAnchorIntersecter.intersect(loopLists, handler,
                chromToHiResAnchorSmoothCounts, highResolution, window);

        mergedDeDuped.exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);
    }

    private static Map<String, float[]> smoothCounts(Map<String, int[]> chromToHiResAnchorCounts) {
        Map<String, float[]> chromToHiResAnchorSmoothCounts = new HashMap<>();
        for (String key : chromToHiResAnchorCounts.keySet()) {
            int[] counts = chromToHiResAnchorCounts.get(key);
            chromToHiResAnchorSmoothCounts.put(key, Convolution1DTools.smooth5(counts));
        }
        return chromToHiResAnchorSmoothCounts;
    }

}
