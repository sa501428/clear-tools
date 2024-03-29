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
// [-filter]
public static String usage = "unwrap2[-legacy] [-r resolution] <genomeID> <loops.bedpe> <output_stem_>\n" +
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

        if (command.contains("2")) {
            Feature2DList invertedList = unwrap2(loopList);
            invertedList.exportFeatureList(new File(outFile + "unwrapped2.local.bedpe"), false, Feature2DList.ListFormat.NA);
        } else {
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
    }

    private static Feature2DList unwrap2(Feature2DList loopList) {
        Feature2DList unwrapped = new Feature2DList();

        loopList.processLists((s, list) -> {
            List<Feature2D> invLocalList = new ArrayList<>(list.size() / 2);

            for (Feature2D feature2D : list) {
                Feature2D inv = unwrapLocal(feature2D);
                if (inv != null) {
                    if (localXYNotNearDiagonal(inv, 5000)) {
                        invLocalList.add(inv);
                    }
                }
            }
            unwrapped.addByKey(s, invLocalList);
        });
        return unwrapped;
    }

    private static Feature2DList[] initializeF2DArray(int n) {
        Feature2DList[] array = new Feature2DList[n];
        for (int k = 0; k < array.length; k++) {
            array[k] = new Feature2DList();
        }
        return array;
    }

    private static Feature2DList[] unwrap(Feature2DList loopList, boolean doFilter, boolean useLegacy) {
        Feature2DList[] unwrapped = initializeF2DArray(numLists);

        loopList.processLists((s, list) -> {
            List<Feature2D> invAnchorsList = new ArrayList<>(list.size() / 2);
            List<Feature2D> invBestAnchorsList = new ArrayList<>(list.size() / 2);
            List<Feature2D> invLocalList = new ArrayList<>(list.size() / 2);
            List<Feature2D> idealList = new ArrayList<>(list.size() / 3);
            List<Feature2D> badAnchorList = new ArrayList<>(list.size() / 3);
            List<Feature2D> badLocalList = new ArrayList<>(list.size() / 3);

            for (Feature2D feature2D : list) {
                Feature2D inv = unwrap(feature2D, useLegacy, true);
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

    private static boolean localXYNotNearDiagonal(Feature2D loop, int minDist) {
        long x = Long.parseLong(loop.getAttribute("localX"));
        long y = Long.parseLong(loop.getAttribute("localY"));
        int dist = (int) Math.min(Math.abs(x - y),
                Math.abs(loop.getMidPt1() - loop.getMidPt2()));
        return dist > minDist;
    }

    private static Feature2D unwrapLocal(Feature2D feature2D) {
        try {
            long start1 = Long.parseLong(feature2D.getAttribute("localX"));
            long start2 = Long.parseLong(feature2D.getAttribute("localY"));
            long end1 = start1 + resolution;
            long end2 = start2 + resolution;
            Map<String, String> attrs = new HashMap<>(feature2D.getAttributes());
            return new Feature2D(Feature2D.FeatureType.PEAK, feature2D.getChr1(), start1, end1,
                    feature2D.getChr2(), start2, end2, Color.BLUE, attrs);
        } catch (Exception ignored) {
        }
        return null;
    }

    private static int getDistanceBetweenAnchorsAndLocal(Feature2D loop) {
        long mid1 = loop.getMidPt1();
        long mid2 = loop.getMidPt2();

        long x = Long.parseLong(loop.getAttribute("localX"));
        long y = Long.parseLong(loop.getAttribute("localY"));

        return (int) (Math.abs(mid1 - x) + Math.abs(mid2 - y));
    }

    public static Feature2D unwrap(Feature2D feature2D, boolean useLegacy, boolean saveOriginalBounds) {
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
            if (saveOriginalBounds) {
                attrs.put("original_start_1", "" + feature2D.getStart1());
                attrs.put("original_start_2", "" + feature2D.getStart2());
                attrs.put("original_end_1", "" + feature2D.getEnd1());
                attrs.put("original_end_2", "" + feature2D.getEnd2());
            }
            return new Feature2D(Feature2D.FeatureType.PEAK, feature2D.getChr1(), start1, end1,
                    feature2D.getChr2(), start2, end2, Color.BLUE, attrs);
        } catch (Exception ignored) {
            return feature2D;
        }
    }
}
