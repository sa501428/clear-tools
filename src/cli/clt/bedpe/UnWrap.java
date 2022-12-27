package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.clean.LoopTools;
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

    public static String usage = "unwrap[-filter][-legacy] [-r resolution] <genomeID> <loops.bedpe> <output_stem_>\n" +
            "\t\tunwrap localizer output to proper hi-res inverted bounds list\n" +
            "\t\tlegacy mode will use the old method of unwrapping the combined anchors\n";
    private static final int numLists = 6;
    private static int resolution = 10;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        resolution = parser.getResolutionOption(10);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);

        Feature2DList loopList = Feature2DParser.loadFeatures(args[2], handler, true, null, false);
        String outFile = args[3];
        boolean doFilter = command.contains("filter");
        boolean useLegacy = command.contains("legacy");

        Feature2DList[] invertedLists = unwrap(loopList, doFilter, useLegacy);
        invertedLists[0].exportFeatureList(new File(outFile + "unwrapped.anchors.bedpe"), false, Feature2DList.ListFormat.NA);
        invertedLists[1].exportFeatureList(new File(outFile + "unwrapped.best.anchors.bedpe"), false, Feature2DList.ListFormat.NA);
        invertedLists[2].exportFeatureList(new File(outFile + "unwrapped.local.bedpe"), false, Feature2DList.ListFormat.NA);

        if (doFilter) {
            invertedLists[3].exportFeatureList(new File(outFile + "optimal.bedpe"), false, Feature2DList.ListFormat.NA);
            invertedLists[4].exportFeatureList(new File(outFile + "suboptimal.anchor.bedpe"), false, Feature2DList.ListFormat.NA);
            invertedLists[5].exportFeatureList(new File(outFile + "suboptimal.local.bedpe"), false, Feature2DList.ListFormat.NA);
        }
    }

    private static Feature2DList[] unwrap(Feature2DList loopList, boolean doFilter, boolean useLegacy) {
        Feature2DList[] unwrapped = new Feature2DList[numLists];
        for (int k = 0; k < unwrapped.length; k++) {
            unwrapped[k] = new Feature2DList();
        }

        loopList.processLists((s, list) -> {
            List<Feature2D> invAnchorsList = new ArrayList<>(list.size() / 2);
            List<Feature2D> invBestAnchorsList = new ArrayList<>(list.size() / 2);
            List<Feature2D> invLocalList = new ArrayList<>(list.size() / 2);
            List<Feature2D> idealList = new ArrayList<>(list.size() / 3);
            List<Feature2D> badAnchorList = new ArrayList<>(list.size() / 3);
            List<Feature2D> badLocalList = new ArrayList<>(list.size() / 3);

            for (Feature2D feature2D : list) {
                Feature2D inv = unwrap(feature2D, useLegacy);
                if (inv != null) {
                    invAnchorsList.add(inv);
                    invLocalList.add(unwrapLocal(feature2D));
                    if (containsLocalAndNotOnDiagonal(inv)) {
                        invBestAnchorsList.add(inv);
                    }

                    int dist = getDistanceBetweenAnchorsAndLocal(inv);
                    inv.addIntAttribute("local_vs_mid_anchor_offset", dist);
                    if (doFilter) {
                        if (dist < 10) {
                            idealList.add(inv);
                        } else if (dist > 400) {
                            badAnchorList.add(inv);
                            badLocalList.add(unwrapLocal(inv));
                        }
                    }
                }
            }
            unwrapped[0].addByKey(s, invAnchorsList);
            unwrapped[1].addByKey(s, invBestAnchorsList);
            unwrapped[2].addByKey(s, invLocalList);
            if (doFilter) {
                unwrapped[3].addByKey(s, idealList);
                unwrapped[4].addByKey(s, badAnchorList);
                unwrapped[5].addByKey(s, badLocalList);
            }
        });
        return unwrapped;
    }

    private static boolean containsLocalAndNotOnDiagonal(Feature2D feature2D) {
        if (LoopTools.dist(feature2D) > 10000) {
            long x = Long.parseLong(feature2D.getAttribute("localX"));
            long y = Long.parseLong(feature2D.getAttribute("localY"));
            return feature2D.getStart1() <= x && x < feature2D.getEnd1() &&
                    feature2D.getStart2() <= y && y < feature2D.getEnd2();
        }
        return false;
    }

    private static Feature2D unwrapLocal(Feature2D feature2D) {
        long start1 = Long.parseLong(feature2D.getAttribute("localX"));
        long start2 = Long.parseLong(feature2D.getAttribute("localY"));
        long end1 = start1 + resolution;
        long end2 = start2 + resolution;
        Map<String, String> attrs = new HashMap<>(feature2D.getAttributes());
        return new Feature2D(Feature2D.FeatureType.PEAK, feature2D.getChr1(), start1, end1,
                feature2D.getChr2(), start2, end2, Color.BLUE, attrs);
    }

    private static int getDistanceBetweenAnchorsAndLocal(Feature2D loop) {
        long mid1 = loop.getMidPt1();
        long mid2 = loop.getMidPt2();

        long x = Long.parseLong(loop.getAttribute("localX"));
        long y = Long.parseLong(loop.getAttribute("localY"));

        return (int) (Math.abs(mid1 - x) + Math.abs(mid2 - y));
    }

    private static Feature2D unwrap(Feature2D feature2D, boolean useLegacy) {
        try {
            long start1, start2, end1, end2;
            if (useLegacy) {
                start1 = Long.parseLong(feature2D.getAttribute("highRes_start_1"));
                start2 = Long.parseLong(feature2D.getAttribute("highRes_start_2"));
                end1 = Long.parseLong(feature2D.getAttribute("highRes_end_1"));
                end2 = Long.parseLong(feature2D.getAttribute("highRes_end_2"));
            } else {
                start1 = Long.parseLong(feature2D.getAttribute("upstream_start_1"));
                end1 = Long.parseLong(feature2D.getAttribute("upstream_end_1"));
                start2 = Long.parseLong(feature2D.getAttribute("downstream_start_2"));
                end2 = Long.parseLong(feature2D.getAttribute("downstream_end_2"));
            }
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
