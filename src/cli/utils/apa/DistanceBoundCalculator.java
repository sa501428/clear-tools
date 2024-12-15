package cli.utils.apa;

import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.util.Collection;
import java.util.List;

public class DistanceBoundCalculator {
    int minDist = Integer.MAX_VALUE;
    int maxDist = 0;
    boolean isIntra;

    private DistanceBoundCalculator(boolean isIntra) {
        this.isIntra = isIntra;
    }

    public static DistanceBoundCalculator fromNestedList(List<List<Feature2D>> loops, int window, int resolution, boolean isIntra) {
        DistanceBoundCalculator calculator = new DistanceBoundCalculator(isIntra);
        if (isIntra) {
            for (List<Feature2D> loop : loops) {
                for (Feature2D f : loop) {
                    int dist = calculator.getDist(f, resolution);
                    calculator.updateDist(dist);
                }
            }
            calculator.adjustBounds(window);
        }
        return calculator;
    }

    public static DistanceBoundCalculator fromFlatList(Collection<Feature2D> loops, int window, int resolution, boolean isIntra) {
        DistanceBoundCalculator calculator = new DistanceBoundCalculator(isIntra);
        if (isIntra) {
            for (Feature2D f : loops) {
                int dist = calculator.getDist(f, resolution);
                calculator.updateDist(dist);
            }
            calculator.adjustBounds(window);
        }
        return calculator;
    }

    private void updateDist(int dist) {
        if (dist < minDist) {
            minDist = dist;
        }
        if (dist > maxDist) {
            maxDist = dist;
        }
    }

    private void adjustBounds(int window) {
        minDist -= 4 * window;
        maxDist += 4 * window;
    }

    private int getDist(Feature2D f, int resolution) {
        return (int) (Math.abs(f.getMidPt1() - f.getMidPt2()) / resolution);
    }

    public boolean inDistanceRange(ContactRecord cr) {
        if (isIntra) {
            int dist = ExpectedUtils.getDist(cr);
            return dist >= minDist && dist <= maxDist;
        } else {
            return true;
        }
    }
}
