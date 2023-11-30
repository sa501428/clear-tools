package cli.utils.clique;

import cli.clt.anchor.AnchorFix;

import java.util.ArrayList;
import java.util.List;

public class Node95 {
    private final List<Long> genomePositions = new ArrayList<>();
    private final List<Long> weakGenomePositions = new ArrayList<>();
    long mu = -1;
    double sigma = -1;

    public Node95(List<Long> longs) {
        genomePositions.addAll(longs);
    }

    public Node95(long x, boolean isSingleton) {
        weakGenomePositions.add(x);
    }

    public static List<Node95> convert(List<List<Long>> clustersToKeep) {
        List<Node95> nodes = new ArrayList<>();
        for (List<Long> cluster : clustersToKeep) {
            nodes.add(new Node95(cluster));
        }
        return nodes;
    }

    public void addWeak(long val) {
        weakGenomePositions.add(val);
    }

    public long[] getBounds95() {
        if (mu < 0 || sigma < 0) {
            calculate95();
        }
        return new long[]{(long) (mu - 1.96 * sigma), (long) (mu + 1.96 * sigma)};
    }

    private double getStdDev(List<Long> genomePositions, long mu) {
        double sum = 0;
        for (long val : genomePositions) {
            sum += ((val - mu) * (val - mu));
        }
        if (sum > 0) {
            return Math.sqrt(sum / genomePositions.size());
        }
        return 0;
    }

    public long getMean(List<Long> genomePositions) {
        long sum = 0;
        for (long val : genomePositions) {
            sum += val;
        }
        return sum / genomePositions.size();
    }

    public long getMu() {
        if (mu == -1) {
            if (genomePositions.size() > 0) {
                mu = getMean(genomePositions);
            } else {
                return getMean(weakGenomePositions);
            }
        }
        return mu;
    }

    public void calculate95() {
        if (genomePositions.size() > 1) {
            mu = getMean(genomePositions);
            sigma = getStdDev(genomePositions, mu);
        } else if (weakGenomePositions.size() > 0) {
            mu = getMean(weakGenomePositions);
            if (weakGenomePositions.size() > 1) {
                sigma = getStdDev(weakGenomePositions, mu);
            } else {
                sigma = AnchorFix.MAX_DIST;
            }
        } else {
            System.err.println("Error: no positions in node");
        }
    }

    public int getWeight() {
        if (genomePositions.size() > 0) {
            return genomePositions.size();
        } else {
            return -weakGenomePositions.size();
        }
    }

    public List<Long> getPositions() {
        return genomePositions;
    }

    public List<Long> getWeakPositions() {
        return weakGenomePositions;
    }
}
