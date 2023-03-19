package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.flags.Anchor;
import cli.utils.flags.AnchorWithScore;
import cli.utils.general.BedFileParser;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FilterBedpeByAnchorAPA {

    public static String usage = "filter-by-anchor-apa [-r resolution] <genomeID> <input.bedpe> <anchors.bed>\n" +
            "\t\tdefault behavior will filter loops by their APA score";
    private static int resolution = 1000;

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(57, usage);
        }
        resolution = parser.getResolutionOption(resolution);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        String infile = args[2];
        Feature2DList loopList = Feature2DParser.loadFeatures(infile, handler, true,
                null, false);
        GenomeWide1DList<Anchor> allAnchors = BedFileParser.loadFromBEDFile(handler, args[3], -1,
                false, true);
        filterByAnchors(loopList, handler, allAnchors, infile);
    }

    private static void filterByAnchors(Feature2DList loopList, ChromosomeHandler handler, GenomeWide1DList<Anchor> allAnchors, String outStem) {
        Feature2DList outputGood = new Feature2DList();
        Feature2DList outputBad = new Feature2DList();
        Feature2DList outputIndeterminate = new Feature2DList();

        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            if (true) System.out.println("Processing " + chrom.getName());
            List<Feature2D> loops = loopList.get(chrom.getIndex(), chrom.getIndex());
            if (loops.size() > 0) {
                List<Anchor> anchors = allAnchors.getFeatures("" + chrom.getIndex());
                if (anchors.size() > 0) {
                    List<Feature2D> goodLoops = new ArrayList<>();
                    List<Feature2D> badLoops = new ArrayList<>();
                    List<Feature2D> indeterminateLoops = new ArrayList<>();
                    filterByLoopAnchors(loops, anchors, goodLoops, badLoops, indeterminateLoops);

                    if (goodLoops.size() > 0) {
                        outputGood.addByKey(Feature2DList.getKey(chrom, chrom), goodLoops);
                    }
                    if (badLoops.size() > 0) {
                        outputBad.addByKey(Feature2DList.getKey(chrom, chrom), badLoops);
                    }
                    if (indeterminateLoops.size() > 0) {
                        outputIndeterminate.addByKey(Feature2DList.getKey(chrom, chrom), indeterminateLoops);
                    }
                }
            }
        }
        outputGood.exportFeatureList(new File(outStem + "_anchorAPA_pass.bed"), false, Feature2DList.ListFormat.NA);
        outputBad.exportFeatureList(new File(outStem + "_anchorAPA_fail.bed"), false, Feature2DList.ListFormat.NA);
        outputIndeterminate.exportFeatureList(new File(outStem + "_anchorAPA_indeterminate.bed"), false, Feature2DList.ListFormat.NA);
    }

    private static void filterByLoopAnchors(List<Feature2D> loops, List<Anchor> anchors, List<Feature2D> goodLoops,
                                            List<Feature2D> badLoops, List<Feature2D> indeterminateLoops) {

        Map<Integer, Anchor> forwardMap = generateAnchorMap(anchors, resolution, true);
        Map<Integer, Anchor> reverseMap = generateAnchorMap(anchors, resolution, false);
        int count = 0;

        // iterate and check for each loop whether the corresponding forward and reverse anchors are available
        for (Feature2D loop : loops) {
            int upStreamBin = (int) (loop.getMidPt1() / resolution);
            int downStreamBin = (int) (loop.getMidPt2() / resolution);

            if (forwardMap.containsKey(upStreamBin) && reverseMap.containsKey(downStreamBin)) {
                AnchorWithScore forwardAnchor = (AnchorWithScore) forwardMap.get(upStreamBin);
                AnchorWithScore reverseAnchor = (AnchorWithScore) reverseMap.get(downStreamBin);

                if (forwardAnchor.getScore() > 1.1 && reverseAnchor.getScore() > 1.1) {
                    goodLoops.add(loop);
                } else if (forwardAnchor.getScore() < 1 || reverseAnchor.getScore() < 1) {
                    badLoops.add(loop);
                } else {
                    indeterminateLoops.add(loop);
                }
            } else {
                System.err.println("Loop anchors not found - something went wrong... (" + (count++) + ")");
            }
        }

    }

    private static Map<Integer, Anchor> generateAnchorMap(List<Anchor> anchors, int resolution, boolean justForwards) {
        Map<Integer, Anchor> anchorMap = new HashMap<>();
        for (Anchor anchor : anchors) {
            AnchorWithScore anchor2 = (AnchorWithScore) anchor;
            if (anchor2.getName().toLowerCase().contains("forward") && justForwards) {
                populateAnchorMap(anchorMap, anchor, resolution);
            } else if (anchor2.getName().toLowerCase().contains("reverse") && !justForwards) {
                populateAnchorMap(anchorMap, anchor, resolution);
            }
        }
        return anchorMap;
    }

    private static void populateAnchorMap(Map<Integer, Anchor> anchorMap, Anchor anchor, int resolution) {
        int width = anchor.getWidth();
        long gx1 = anchor.getStart();
        int bin = (int) ((gx1 + width) / resolution);
        anchorMap.put(bin, anchor);
    }
}
