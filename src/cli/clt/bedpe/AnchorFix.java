package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.clique.Node95;
import cli.utils.clique.SimpleClustering;
import cli.utils.flags.Anchor;
import cli.utils.general.BedTools;
import cli.utils.loops.AnchorTools;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.*;

public class AnchorFix {

    public static String usage = "anchor-fix[-clean] [-r resolution] <genomeID> <input.bedpe> <output.stem>\n" +
            "\t\tdefault behavior will fix the hi-res shared anchors for loops\n" +
            "\t\tclean avoids saving old attributes";
    private static int resolution = 100;
    public static int MAX_DIST = 250;

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
                List<Feature2D> newLoops = recoverLoops(loops);
                output.addByKey(Feature2DList.getKey(chrom, chrom), newLoops);
            }
        }

        output.exportFeatureList(new File(outStem + ".fixed.anchors.bedpe"), false, Feature2DList.ListFormat.NA);
        exportAnchors(output, outStem);
    }

    private static void exportAnchors(Feature2DList output, String outStem) {
        Set<Anchor> disoriented = new HashSet<>();
        Set<Anchor> upstream = new HashSet<>();
        Set<Anchor> downstream = new HashSet<>();
        output.processLists((s, list) -> {
            for (Feature2D f : list) {
                disoriented.add(AnchorTools.getAnchor(f, "highRes_start_1", "highRes_end_1"));
                disoriented.add(AnchorTools.getAnchor(f, "highRes_start_2", "highRes_end_2"));
                upstream.add(AnchorTools.getAnchor(f, "upstream_start_1", "upstream_end_1"));
                downstream.add(AnchorTools.getAnchor(f, "downstream_start_2", "downstream_end_2"));
            }
        });

        BedTools.exportBedFile(new File(outStem + ".disoriented.anchors.bed"), disoriented);
        BedTools.exportBedFile(new File(outStem + ".upstream.anchors.bed"), upstream);
        BedTools.exportBedFile(new File(outStem + ".downstream.anchors.bed"), downstream);
    }

    private static List<Feature2D> recoverLoops(List<Feature2D> loops) {
        List<Long> upStreamAnchorBins = getAllOfFeature(loops, "localX");
        List<Long> downStreamAnchorBins = getAllOfFeature(loops, "localY");
        List<Long> allAnchorBins = new ArrayList<>(2 * loops.size());
        allAnchorBins.addAll(upStreamAnchorBins);
        allAnchorBins.addAll(downStreamAnchorBins);

        List<Node95> upStreamNodes = getNodes(upStreamAnchorBins, resolution);
        List<Node95> downStreamNodes = getNodes(downStreamAnchorBins, resolution);
        List<Node95> allNodes = getNodes(allAnchorBins, resolution);

        return fixedList(loops, upStreamNodes, downStreamNodes, allNodes);
    }

    private static List<Feature2D> fixedList(List<Feature2D> loops, List<Node95> upStreamNodes,
                                             List<Node95> downStreamNodes, List<Node95> allNodes) {
        Map<Long, Node95> upStreamBinToNode = buildPositionToNodeMapping(upStreamNodes);
        Map<Long, Node95> downStreamBinToNode = buildPositionToNodeMapping(downStreamNodes);
        Map<Long, Node95> anyBinToNode = buildPositionToNodeMapping(allNodes);

        for (Feature2D loop : loops) {
            long x = Long.parseLong(loop.getAttribute("localX"));
            long y = Long.parseLong(loop.getAttribute("localY"));
            Node95 xNodeUp = upStreamBinToNode.get(x);
            Node95 yNodeDown = downStreamBinToNode.get(y);
            Node95 xNodeAny = anyBinToNode.get(x);
            Node95 yNodeAny = anyBinToNode.get(y);
            if (xNodeUp == null || yNodeDown == null || xNodeAny == null || yNodeAny == null) {
                System.err.println("Error: null node :  " + x + "  " + y + "\n " + xNodeUp + "\n" + yNodeDown + "\n" + xNodeAny + "\n" + yNodeAny);
            }
            setAttributes(loop, xNodeUp, yNodeDown, xNodeAny, yNodeAny);
        }
        return loops;
    }

    private static void setAttributes(Feature2D loop, Node95 xNodeUp, Node95 yNodeDown,
                                      Node95 xNodeAny, Node95 yNodeAny) {
        if (xNodeUp == null || yNodeDown == null || xNodeAny == null || yNodeAny == null) {
            //System.err.println("Error: null node\n "+xNodeUp+"\n"+yNodeDown+"\n"+xNodeAny+"\n"+yNodeAny);
        } else {
            setAttributes(loop, xNodeUp, "upstream_start_1", "upstream_end_1");
            setAttributes(loop, yNodeDown, "downstream_start_2", "downstream_end_2");
            setAttributes(loop, xNodeAny, "highRes_start_1", "highRes_end_1");
            setAttributes(loop, yNodeAny, "highRes_start_2", "highRes_end_2");
        }
    }

    private static void setAttributes(Feature2D loop, Node95 node, String startName, String endName) {
        if (node == null) {
            System.err.println("Node is null _ " + startName + " - " + endName);
        }
        long[] bounds = node.getBounds95();
        loop.setAttribute(startName, String.valueOf(bounds[0]));
        loop.setAttribute(endName, String.valueOf(bounds[1]));
    }


    private static List<Long> getAllOfFeature(List<Feature2D> loops, String attribute) {
        List<Long> positions = new ArrayList<>(loops.size());
        for (Feature2D feature : loops) {
            positions.add(Long.parseLong(feature.getAttribute(attribute)));
        }
        return positions;
    }


    public static List<Node95> getNodes(List<Long> genomePositions, int resolution) {
        List<List<Long>> clusters = SimpleClustering.cluster(genomePositions, resolution);
        List<Long> pointsToReassign = new ArrayList<>(clusters.size() / 2);
        List<List<Long>> clustersToKeep = new ArrayList<>(clusters.size() / 2);
        for (List<Long> cluster : clusters) {
            if (cluster.size() > 1) {
                clustersToKeep.add(cluster);
            } else {
                pointsToReassign.add(cluster.get(0));
            }
        }

        List<Node95> cleanedUpClusters = assignToNearestClusterOrMakeSingleton(clustersToKeep, pointsToReassign);
        clustersToKeep.clear();
        pointsToReassign.clear();
        return cleanedUpClusters;
    }

    private static List<Node95> assignToNearestClusterOrMakeSingleton(List<List<Long>> clustersToKeep, List<Long> pointsToReassign) {
        List<Node95> nodes = Node95.convert(clustersToKeep);
        for (Node95 node : nodes) {
            node.calculate95();
        }
        for (long point : pointsToReassign) {
            Node95 nearest = getNearestNode(point, nodes, MAX_DIST);
            if (nearest == null) {
                nodes.add(new Node95(point, true));
            } else {
                nearest.addWeak(point);
            }
        }
        return nodes;
    }

    private static Node95 getNearestNode(long point, List<Node95> nodes, int maxDist) {
        double currDist = Double.MAX_VALUE;
        Node95 nearest = null;
        for (Node95 node : nodes) {
            double dist = Math.abs(node.getMu() - point);
            if (dist < currDist) {
                currDist = dist;
                nearest = node;
            }
        }
        if (currDist < maxDist) {
            return nearest;
        }
        return null;
    }

    public static Map<Long, Node95> buildPositionToNodeMapping(List<Node95> nodes) {
        Map<Long, Node95> mapping = new HashMap<>();
        for (Node95 node : nodes) {
            for (long val : node.getPositions()) {
                mapping.put(val, node);
            }
            for (long val : node.getWeakPositions()) {
                mapping.put(val, node);
            }
        }
        return mapping;
    }
}
