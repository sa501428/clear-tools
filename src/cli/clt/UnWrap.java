package cli.clt;

import cli.Main;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class UnWrap {

    public static String usage = "unwrap <genomeID> <loops.bedpe> <output.bedpe>\n" +
            "\t\tunwrap localizer output to proper hi-res inverted bounds list";

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);

        Feature2DList loopList = Feature2DParser.loadFeatures(args[2], handler, true, null, false);
        String outFile = args[3];

        Feature2DList invertedList = unwrap(loopList);
        invertedList.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static Feature2DList unwrap(Feature2DList loopList) {
        Feature2DList unwrapped = new Feature2DList();
        loopList.processLists((s, list) -> {
            List<Feature2D> invList = new ArrayList<>();
            for (Feature2D feature2D : list) {
                Feature2D inv = unwrap(feature2D);
                if (inv != null) {
                    invList.add(inv);
                }
            }
            unwrapped.addByKey(s, invList);
        });
        return unwrapped;
    }

    private static Feature2D unwrap(Feature2D feature2D) {
        try {
            long start1 = Long.parseLong(feature2D.getAttribute("highRes_start_1"));
            long start2 = Long.parseLong(feature2D.getAttribute("highRes_start_2"));
            long end1 = Long.parseLong(feature2D.getAttribute("highRes_end_1"));
            long end2 = Long.parseLong(feature2D.getAttribute("highRes_end_2"));
            Map<String, String> attrs = new HashMap<>(feature2D.getAttributes());
            attrs.put("original_start_1", "" + feature2D.getStart1());
            attrs.put("original_start_2", "" + feature2D.getStart2());
            attrs.put("original_end_1", "" + feature2D.getEnd1());
            attrs.put("original_end_2", "" + feature2D.getEnd2());
            return new Feature2D(Feature2D.FeatureType.PEAK, feature2D.getChr1(), start1, end1,
                    feature2D.getChr2(), start2, end2, Color.BLUE, attrs);
        } catch (Exception ignored) {
        }
        return null;
    }
}
