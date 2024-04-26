package cli.utils.apa;

import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.util.List;

public class DistanceBoundCalculator {
    int minDist = Integer.MAX_VALUE;
    int maxDist = 0;

    public DistanceBoundCalculator(List<List<Feature2D>> loops, int window, int resolution) {

        for (List<Feature2D> loop : loops) {
            for (Feature2D f : loop) {
                int dist = getDist(f, resolution);
                if (dist < minDist) {
                    minDist = dist;
                }
                if (dist > maxDist) {
                    maxDist = dist;
                }
            }
        }

        minDist -= 4 * window;
        maxDist += 4 * window;
    }

    private int getDist(Feature2D f, int resolution) {
        return (int) (Math.abs(f.getMidPt1() - f.getMidPt2()) / resolution);
    }

    public boolean inDistanceRange(ContactRecord cr) {
        int dist = ExpectedUtils.getDist(cr);
        return dist >= minDist && dist <= maxDist;
    }
}
