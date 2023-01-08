package cli.utils.clique;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class Unused {
    public static void run(int[] counts, AtomicInteger nodeCount) {
        while (hasNonZero(counts)) {
            int index = getMaxIndex(counts);
            Node newNode = new Node(nodeCount.getAndIncrement(), index, new ArrayList<>());
            counts[index] = 0;
            expandInDirection(index + 1, 1, counts, newNode);
            expandInDirection(index - 1, -1, counts, newNode);
        }
    }

    private static void expandInDirection(int start, int increment, int[] counts, Node newNode) {
        int i = start;
        while (true) {
            if (i < 0 || i > counts.length - 1) break;
            if (counts[i] > 0) {
                newNode.add(i, new ArrayList<>());
                counts[i] = 0;
            } else {
                break;
            }
            i = i + increment;
        }
    }

    private static int getMaxIndex(int[] counts) {
        int maxIndex = 0;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > counts[maxIndex]) {
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    private static boolean hasNonZero(int[] counts) {
        for (int val : counts) {
            if (val > 0) return true;
        }
        return false;
    }


    public static List<Node> getNodes(List<Long> genomePositions, AtomicInteger nodeCount, int resolution) {
        Collections.sort(genomePositions);
        int maxBin = (int) ((genomePositions.get(genomePositions.size() - 1) / resolution) + 1);

        int[] counts = getCounts(genomePositions, resolution);

        List<Node> nodes = new ArrayList<>(genomePositions.size() / 2);

        Node current = null;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > 0) {
                if (current == null) {
                    current = new Node(nodeCount.getAndIncrement(), i, new ArrayList<>());
                    nodes.add(current);
                } else {
                    current.add(i, new ArrayList<>());
                }
            } else {
                current = null;
            }
        }

        return nodes;
    }


    private static int[] getCounts(List<Long> genomePositions, int resolution) {
        int maxBin = (int) ((genomePositions.get(genomePositions.size() - 1) / resolution) + 1);
        int[] counts = new int[maxBin];
        for (long val : genomePositions) {
            counts[(int) (val / resolution)]++;
        }
        return counts;
    }

}
