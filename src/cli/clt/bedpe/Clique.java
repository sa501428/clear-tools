package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.clique.MatrixUtils;
import cli.utils.clique.NetworkBuilder;
import cli.utils.clique.Node;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class Clique {

    // [-rescue]
    public static String usage = "clique[-clean][-rescue] [-r resolution] <genomeID> <input.bedpe> <output.stem>\n" +
            "\t\tdefault behavior finds the cliques using midpoints of the anchors at the resolution specified\n" +
            "\t\trescue will predict loops that were potentially missed\n" +
            "\t\tclean avoids saving old attributes";
    private static int resolution = 200;

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }
        resolution = parser.getResolutionOption(200);
        String genomeID = args[1];
        String inFile = args[2];
        String outStem = args[3];
        if (command.contains("rescue")) {
            rescueLoops(inFile, genomeID, outStem, command.contains("clean"));
            System.out.println("clique rescue complete");
        } else {
            // todo
        }
    }

    private static void rescueLoops(String inputBedpe, String genomeID, String outStem, boolean noAttributes) {
        Feature2DList output = new Feature2DList();
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpe, handler, !noAttributes,
                null, false);
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
            List<Feature2D> loops = loopList.get(chrom.getIndex(), chrom.getIndex());
            if (loops.size() > 0) {
                if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
                List<Feature2D> newLoops = recoverLoops(loops);
                output.addByKey(Feature2DList.getKey(chrom, chrom), newLoops);
            }
        }

        output.exportFeatureList(new File(outStem + ".rescue.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> recoverLoops(List<Feature2D> initialLoops) {

        List<Long> upStreamAnchorBins = new ArrayList<>(initialLoops.size());
        List<Long> downStreamAnchorBins = new ArrayList<>(initialLoops.size());
        for (Feature2D feature : initialLoops) {
            upStreamAnchorBins.add((feature.getMidPt1()));
            downStreamAnchorBins.add((feature.getMidPt2()));
        }

        AtomicInteger nodeCount = new AtomicInteger(0);
        List<Node> upStreamNodes = NetworkBuilder.getNodes(upStreamAnchorBins, nodeCount, resolution);
        List<Node> downStreamNodes = NetworkBuilder.getNodes(downStreamAnchorBins, nodeCount, resolution);
        int maxN = nodeCount.get();

        Map<Integer, Node> upStreamBinToNode = NetworkBuilder.buildIndexToNodeMapping(upStreamNodes);
        Map<Integer, Node> downStreamBinToNode = NetworkBuilder.buildIndexToNodeMapping(downStreamNodes);

        float[][] adjacencyMatrix = NetworkBuilder.buildAdjacencyMatrix(maxN, initialLoops,
                upStreamBinToNode, downStreamBinToNode, resolution);
        float[][] a3 = MatrixUtils.cube(adjacencyMatrix);
        MatrixUtils.setEverythingBelowDiagonalToZero(a3);
        adjacencyMatrix = null;

        Map<Integer, Node> idToNode = new HashMap<>();
        for (Node node : upStreamNodes) {
            idToNode.put(node.getId(), node);
        }
        for (Node node : downStreamNodes) {
            idToNode.put(node.getId(), node);
        }

        return retrieveAllLoopsPerAdjMatrix(a3, idToNode, initialLoops.get(0).getChr1(), initialLoops.size());
    }

    private static List<Feature2D> retrieveAllLoopsPerAdjMatrix(float[][] adjMatrix, Map<Integer, Node> idToNode,
                                                                String chrom, int numLoops) {
        List<Feature2D> newLoops = new ArrayList<>(2 * numLoops);
        for (int i = 0; i < adjMatrix.length; i++) {
            for (int j = i + 1; j < adjMatrix[0].length; j++) {
                if (adjMatrix[i][j] > 0) {
                    Node node1 = idToNode.get(i);
                    Node node2 = idToNode.get(j);
                    if (node1.getMinPosition() < node2.getMinPosition()) {
                        newLoops.add(new Feature2D(Feature2D.FeatureType.PEAK,
                                chrom, node1.getMinPosition() - resolution, node1.getMaxPosition() + resolution,
                                chrom, node2.getMinPosition() - resolution, node2.getMaxPosition() + resolution,
                                Color.BLACK, new HashMap<>()));
                    } else {
                        newLoops.add(new Feature2D(Feature2D.FeatureType.PEAK,
                                chrom, node2.getMinPosition() - resolution, node2.getMaxPosition() + resolution,
                                chrom, node1.getMinPosition() - resolution, node1.getMaxPosition() + resolution,
                                Color.BLACK, new HashMap<>()));
                    }
                }
            }
        }
        return newLoops;
    }
}
