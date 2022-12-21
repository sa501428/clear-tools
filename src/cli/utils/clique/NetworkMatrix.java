package cli.utils.clique;

import javastraw.feature2D.Feature2D;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class NetworkMatrix {
    private final float[][] adjacencyMatrix;
    private final List<Node> upStreamNodes;
    private final List<Node> downStreamNodes;
    private final Map<Integer, Node> idToNode;

    public NetworkMatrix(List<Feature2D> loops, int resolution) {
        List<Long> upStreamAnchorBins = new ArrayList<>(loops.size());
        List<Long> downStreamAnchorBins = new ArrayList<>(loops.size());
        for (Feature2D feature : loops) {
            upStreamAnchorBins.add((feature.getMidPt1()));
            downStreamAnchorBins.add((feature.getMidPt2()));
        }

        AtomicInteger nodeCount = new AtomicInteger(0);
        this.upStreamNodes = NetworkBuilder.getNodes(upStreamAnchorBins, nodeCount, resolution);
        this.downStreamNodes = NetworkBuilder.getNodes(downStreamAnchorBins, nodeCount, resolution);
        int maxN = nodeCount.get();

        Map<Integer, Node> upStreamBinToNode = NetworkBuilder.buildIndexToNodeMapping(upStreamNodes);
        Map<Integer, Node> downStreamBinToNode = NetworkBuilder.buildIndexToNodeMapping(downStreamNodes);

        this.adjacencyMatrix = NetworkBuilder.buildAdjacencyMatrix(maxN, loops,
                upStreamBinToNode, downStreamBinToNode, resolution);

        this.idToNode = buildIDToNodeMapping(upStreamNodes, downStreamNodes);

    }

    private static Map<Integer, Node> buildIDToNodeMapping(List<Node> upStreamNodes, List<Node> downStreamNodes) {
        Map<Integer, Node> idToNode = new HashMap<>();
        for (Node node : upStreamNodes) {
            idToNode.put(node.getId(), node);
        }
        for (Node node : downStreamNodes) {
            idToNode.put(node.getId(), node);
        }
        return idToNode;
    }

    public float[][] getAdjMatrix() {
        return adjacencyMatrix;
    }

    public Map<Integer, Node> getIDToNodeMapping() {
        return idToNode;
    }
}
