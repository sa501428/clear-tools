package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.clique.ConnectedComponents;
import cli.utils.clique.MatrixUtils;
import cli.utils.clique.NetworkMatrix;
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
            splitFilesIntoCliques(inFile, genomeID, outStem, command.contains("clean"));
        }
    }

    private static void splitFilesIntoCliques(String inputBedpe, String genomeID, String outStem, boolean noAttributes) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpe, handler, !noAttributes,
                null, false);
        int cliqueCounter = 0;
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
            List<Feature2D> loops = loopList.get(chrom.getIndex(), chrom.getIndex());
            if (loops.size() > 0) {
                if (Main.printVerboseComments) System.out.println("Processing " + chrom.getName());
                List<List<Feature2D>> cliques = getCliques(loops);
                for (List<Feature2D> clique : cliques) {
                    if (clique.size() > 0) {
                        Feature2DList output = new Feature2DList();
                        output.addByKey(Feature2DList.getKey(chrom, chrom), clique);
                        output.exportFeatureList(new File(outStem + ".clique." + cliqueCounter + ".bedpe"),
                                false, Feature2DList.ListFormat.NA);
                        cliqueCounter++;
                    }
                }
            }
        }
    }

    private static List<List<Feature2D>> getCliques(List<Feature2D> loops) {
        NetworkMatrix networkMatrix = new NetworkMatrix(loops, resolution);
        return retrieveAllCliquesPerAdjMatrix(networkMatrix.getAdjMatrix(),
                networkMatrix.getIDToNodeMapping(), loops.get(0).getChr1());
    }

    private static List<List<Feature2D>> retrieveAllCliquesPerAdjMatrix(float[][] adjMatrix, Map<Integer, Node> idToNodeMapping, String chr1) {
        MatrixUtils.setEverythingBelowDiagonalToZero(adjMatrix);
        List<List<Feature2D>> cliques = new ArrayList<>();
        for (int i = 0; i < adjMatrix.length; i++) {
            if (hasNonZeroEntries(adjMatrix[i])) {
                List<Feature2D> clique = ConnectedComponents.processClique(adjMatrix, idToNodeMapping, i, chr1, resolution);
                cliques.add(clique);
            }
        }
        return cliques;
    }

    private static boolean hasNonZeroEntries(float[] row) {
        for (float entry : row) {
            if (entry > 0) {
                return true;
            }
        }
        return false;
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
        NetworkMatrix matrix = new NetworkMatrix(initialLoops, resolution);
        float[][] a3 = MatrixUtils.cube(matrix.getAdjMatrix());
        MatrixUtils.setEverythingBelowDiagonalToZero(a3);
        Map<Integer, Node> idToNode = matrix.getIDToNodeMapping();
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
