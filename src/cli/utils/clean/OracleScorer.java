package cli.utils.clean;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class OracleScorer {
    public static Feature2DList filter(Feature2DList loopList, double threshold) {
        loopList.filterLists((chr, feature2DList) -> filterByOracleParams(feature2DList, threshold));
        return loopList;
    }

    private static List<Feature2D> filterByOracleParams(List<Feature2D> loops, double threshold) {
        Set<Feature2D> goodLoops = new HashSet<>();
        for (Feature2D loop : loops) {
            if (passesOracleCriteria(loop, threshold)) {
                goodLoops.add(loop);
            }
        }
        return new ArrayList<>(goodLoops);
    }

    private static boolean passesOracleCriteria(Feature2D loop, double cutoff) {
        float r1k = Float.parseFloat(loop.getAttribute("score_oracle_1000"));
        float r2k = Float.parseFloat(loop.getAttribute("score_oracle_2000"));
        float r5k = Float.parseFloat(loop.getAttribute("score_oracle_5000"));
        float r10k = Float.parseFloat(loop.getAttribute("score_oracle_10000"));

        return (r1k > cutoff) || (r2k > cutoff) ||
                ((r1k > cutoff) && (r2k > cutoff)) ||
                ((r2k > cutoff) && (r5k > cutoff)) ||
                ((r1k > cutoff) && (r2k > cutoff) && (r5k > cutoff)) ||
                ((r1k > cutoff) && (r2k > cutoff) && (r10k > cutoff)) ||
                ((r1k > cutoff) && (r5k > cutoff) && (r10k > cutoff)) ||
                ((r2k > cutoff) && (r5k > cutoff) && (r10k > cutoff));
    }
}
