package cli.utils.clique;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class SimpleClustering {
    public static List<List<Long>> cluster(List<Long> genomePositions, int resolution) {
        Collections.sort(genomePositions);
        int maxBin = (int) ((genomePositions.get(genomePositions.size() - 1) / resolution) + 1);

        Map<Integer, List<Long>> counts = NetworkBuilder.buildBinToPositionsMapping(genomePositions, resolution);

        List<List<Long>> nodes = new ArrayList<>(genomePositions.size() / 2);

        List<Long> current = null;
        for (int i = 0; i < maxBin; i++) {
            if (counts.containsKey(i)) {
                if (current == null) {
                    current = counts.get(i);
                    nodes.add(current);
                } else {
                    current.addAll(counts.get(i));
                }
            } else {
                current = null;
            }
        }

        return nodes;
    }
}
