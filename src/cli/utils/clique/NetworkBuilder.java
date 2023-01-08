package cli.utils.clique;

import javastraw.feature2D.Feature2D;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class NetworkBuilder {

    public static List<Node> getNodes(List<Long> genomePositions, AtomicInteger nodeCount, int resolution) {
        Collections.sort(genomePositions);
        int maxBin = (int) ((genomePositions.get(genomePositions.size() - 1) / resolution) + 1);

        Map<Integer, List<Long>> counts = buildBinToPositionsMapping(genomePositions, resolution);

        List<Node> nodes = new ArrayList<>(genomePositions.size() / 2);

        Node current = null;
        for (int i = 0; i < maxBin; i++) {
            if (counts.containsKey(i)) {
                if (current == null) {
                    current = new Node(nodeCount.getAndIncrement(), i, counts.get(i));
                    nodes.add(current);
                } else {
                    current.add(i, counts.get(i));
                }
            } else {
                current = null;
            }
        }

        return nodes;
    }

    public static Map<Integer, List<Long>> buildBinToPositionsMapping(List<Long> genomePositions, int resolution) {
        Map<Integer, List<Long>> mapping = new HashMap<>();
        for (long val : genomePositions) {
            int bin = (int) (val / resolution);
            if (!mapping.containsKey(bin)) {
                mapping.put(bin, new LinkedList<>());
            }
            mapping.get(bin).add(val);
        }
        return mapping;
    }

    public static Map<Integer, Node> buildIndexToNodeMapping(List<Node> nodes) {
        Map<Integer, Node> mapping = new HashMap<>();
        for (Node node : nodes) {
            for (int i = node.getMinIndex(); i < node.getMaxIndex() + 1; i++) {
                mapping.put(i, node);
            }
        }
        return mapping;
    }

    public static float[][] buildAdjacencyMatrix(int maxN, List<Feature2D> list,
                                                 Map<Integer, Node> upStreamBinToNode,
                                                 Map<Integer, Node> downStreamBinToNode,
                                                 int resolution) {
        float[][] matrix = new float[maxN][maxN];
        for (Feature2D feature : list) {
            int upStreamBin = (int) (feature.getMidPt1() / resolution);
            int downStreamBin = (int) (feature.getMidPt2() / resolution);
            if (upStreamBinToNode.containsKey(upStreamBin) && downStreamBinToNode.containsKey(downStreamBin)) {
                int upStreamID = upStreamBinToNode.get(upStreamBin).getId();
                int downStreamID = downStreamBinToNode.get(downStreamBin).getId();
                if (upStreamID != downStreamID) {
                    matrix[upStreamID][downStreamID] = 1;
                    matrix[downStreamID][upStreamID] = 1;
                }
            }
        }
        return matrix;
    }
}
