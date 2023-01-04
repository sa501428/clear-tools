package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.clique.Node95;
import cli.utils.peaks.Point1D;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;

import java.io.File;
import java.util.*;

public class AnchorFix {

    private static final int widthToConnect = 2;
    public static String usage = "anchor-fix[-clean] [-r resolution] <genomeID> <input.bedpe> <output.stem>\n" +
            "\t\tdefault behavior will fix the hi-res shared anchors for loops\n" +
            "\t\tclean avoids saving old attributes";
    private static int resolution = 100;
    public static int MAX_DIST = 250;
    public static int CLUSTER_DIST = 100;

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
                if (true) System.out.println("Processing " + chrom.getName());
                List<Feature2D> newLoops = recoverLoops(loops);
                output.addByKey(Feature2DList.getKey(chrom, chrom), newLoops);
            }
        }

        output.exportFeatureList(new File(outStem + ".fixed.anchors.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> recoverLoops(List<Feature2D> loops) {
        List<Long> upStreamAnchorBins = getAllOfFeature(loops, "localX");
        List<Long> downStreamAnchorBins = getAllOfFeature(loops, "localY");
        List<Long> allAnchorBins = new ArrayList<>(2 * loops.size());
        allAnchorBins.addAll(upStreamAnchorBins);
        allAnchorBins.addAll(downStreamAnchorBins);

        List<Node95> upStreamNodes = getNodes(upStreamAnchorBins);
        List<Node95> downStreamNodes = getNodes(downStreamAnchorBins);
        List<Node95> allNodes = getNodes(allAnchorBins);

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


    public static List<Node95> getNodes(List<Long> genomePositions) {

        /*
        List<Cluster<Point1D>> clusters = dbscan1D(genomePositions);
        List<Point1D> pointsToReassign = new ArrayList<>(clusters.size()/2);
        List<Cluster<Point1D>> clustersToKeep = new ArrayList<>(clusters.size()/2);
        for(Cluster<Point1D> cluster : clusters) {
            if(cluster.getPoints().size() > 1){
                clustersToKeep.add(cluster);
            } else {
                pointsToReassign.add(cluster.getPoints().get(0));
            }
        }

        List<Node95> cleanedUpClusters = assignToNearestClusterOrMakeSingleton(clustersToKeep, pointsToReassign);
        clustersToKeep.clear();
        pointsToReassign.clear();
        return cleanedUpClusters;
        */
        List<Node95> test = new ArrayList<>(genomePositions.size());
        for (Long pos : genomePositions) {
            test.add(new Node95(pos, false));
        }
        return test;
    }

    private static List<Node95> assignToNearestClusterOrMakeSingleton(List<Cluster<Point1D>> clustersToKeep, List<Point1D> pointsToReassign) {
        List<Node95> nodes = Node95.convert(clustersToKeep);
        for (Node95 node : nodes) {
            node.calculate95();
        }
        for (Point1D point : pointsToReassign) {
            Node95 nearest = getNearestNode(point, nodes, MAX_DIST);
            if (nearest == null) {
                nodes.add(new Node95(point.getX(), true));
            } else {
                nearest.addWeak(point.getX());
            }
        }
        return nodes;
    }

    private static Node95 getNearestNode(Point1D point, List<Node95> nodes, int maxDist) {
        double currDist = Double.MAX_VALUE;
        Node95 nearest = null;
        for (Node95 node : nodes) {
            double dist = Math.abs(node.getMu() - point.getX());
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

    private static List<Cluster<Point1D>> dbscan1D(List<Long> points) {
        DBSCANClusterer<Point1D> dbscan = new DBSCANClusterer<>(CLUSTER_DIST, 1);
        return dbscan.cluster(convert(points));
    }

    private static Collection<Point1D> convert(List<Long> points) {
        List<Point1D> converted = new ArrayList<>(points.size());
        for (Long point : points) {
            converted.add(new Point1D(point));
        }
        return converted;
    }

    private static int[] getCounts(Map<Integer, List<Long>> counts, int maxBin) {
        int[] countsTrack = new int[maxBin];
        for (int i = 0; i < maxBin; i++) {
            if (counts.containsKey(i)) {
                countsTrack[i] = counts.get(i).size();
            }
        }
        return countsTrack;
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
