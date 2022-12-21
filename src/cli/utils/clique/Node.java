package cli.utils.clique;

import java.util.ArrayList;
import java.util.List;

public class Node {
    private final int id;
    private final List<Long> genomePositions = new ArrayList<>();
    private int minIndex;
    private int maxIndex;

    public Node(int id, int index, List<Long> longs) {
        this.id = id;
        minIndex = index;
        maxIndex = index;
        genomePositions.addAll(longs);
    }

    void add(int index, List<Long> longs) {
        minIndex = Math.min(minIndex, index);
        maxIndex = Math.max(maxIndex, index);
        genomePositions.addAll(longs);
    }

    public int getId() {
        return id;
    }

    public int getMinIndex() {
        return minIndex;
    }

    public int getMaxIndex() {
        return maxIndex;
    }

    public long getMinPosition() {
        long minVal = genomePositions.get(0);
        for (long val : genomePositions) {
            minVal = Math.min(minVal, val);
        }
        return minVal;
    }

    public long getMaxPosition() {
        long maxVal = genomePositions.get(0);
        for (long val : genomePositions) {
            maxVal = Math.max(maxVal, val);
        }
        return maxVal;
    }
}
