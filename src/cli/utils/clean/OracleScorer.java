package cli.utils.clean;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class OracleScorer {
    public static Feature2DList filter(Feature2DList loopList, boolean beStrict) {
        loopList.filterLists((chr, feature2DList) -> filterByOracleParams(feature2DList, beStrict));
        return loopList;
    }

    private static List<Feature2D> filterByOracleParams(List<Feature2D> loops, boolean beStrict) {
        Set<Feature2D> goodLoops = new HashSet<>();
        for (Feature2D loop : loops) {
            if (passesOracleCriteria(loop, beStrict)) {
                goodLoops.add(loop);
            }
        }
        return new ArrayList<>(goodLoops);
    }

    private static boolean passesOracleCriteria(Feature2D loop, boolean beStrict) {
        float r1k = Float.parseFloat(loop.getAttribute("score_oracle_1000"));
        float r2k = Float.parseFloat(loop.getAttribute("score_oracle_2000"));
        float r5k = Float.parseFloat(loop.getAttribute("score_oracle_5000"));
        float r10k = Float.parseFloat(loop.getAttribute("score_oracle_10000"));

        int dist = LoopTools.dist(loop);
        if (dist > 1000000) {
            if (r5k > 0.85 || r10k > 0.85) {
                if (beStrict) {
                    return r5k > 0.85 && r10k > 0.85;
                } else {
                    return true;
                }
            }
        }

        if (r1k > 0.2 || r2k > 0.2) {
            if (beStrict) {
                return r1k > 0.2 && r2k > 0.2;
            } else {
                return true;
            }
        }
        return false;
    }

    private static boolean defaultCheck4(float r1k, float r2k, float r5k, float r10k, double cutoff, int dist) {
        boolean lowResCheck = false;
        if (dist > 1000000) {
            lowResCheck = ((r5k > cutoff) && (r10k > cutoff));
        }

        return (r1k > cutoff) || (r2k > cutoff);//||
        //((r1k > cutoff) && (r2k > cutoff)) ||
        //((r2k > cutoff) && (r5k > cutoff)) ||
        //lowResCheck ||
        //((r1k > cutoff) && (r2k > cutoff) && (r5k > cutoff)) ||
        //((r1k > cutoff) && (r2k > cutoff) && (r10k > cutoff)) ||
        //((r1k > cutoff) && (r5k > cutoff) && (r10k > cutoff)) ||
        //((r2k > cutoff) && (r5k > cutoff) && (r10k > cutoff));
    }

    private static boolean defaultCheck3(float r1k, float r2k, float r5k, double cutoff) {
        return (r1k > cutoff) || (r2k > cutoff) ||
                ((r1k > cutoff) && (r2k > cutoff)) ||
                ((r2k > cutoff) && (r5k > cutoff)) ||
                ((r1k > cutoff) && (r2k > cutoff) && (r5k > cutoff));
    }

    private static boolean defaultCheck2(float r1k, float r2k, double cutoff) {
        return (r1k > cutoff) || (r2k > cutoff) ||
                ((r1k > cutoff) && (r2k > cutoff));
    }
}
