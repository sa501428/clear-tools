package cli.utils.clean;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DFilter;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;

import java.util.ArrayList;
import java.util.List;

public class LoopTools {

    public static Feature2DList loadNearDiagonalFilteredBedpe(String bedpeFile, ChromosomeHandler handler,
                                                              boolean loadAttributes) {
        return Feature2DParser.loadFeatures(bedpeFile, handler,
                loadAttributes, getNearDiagonalFilter(), false);
    }

    private static Feature2DFilter getNearDiagonalFilter() {
        return (s, initialLoops) -> {
            List<Feature2D> result = new ArrayList<>(initialLoops.size());
            for (Feature2D loop : initialLoops) {
                if (passesMinLoopSize(loop)) {
                    result.add(loop);
                }
            }
            return result;
        };
    }

    public static boolean passesMinLoopSize(Feature2D loop) {
        return dist(loop) / Math.max(getResolutionLoopWasCalledAt(loop), 1000) > 5;
    }

    public static int getResolutionLoopWasCalledAt(Feature2D loop) {
        return (int) Math.max(Math.max(loop.getWidth1(), loop.getWidth2()), 1);
    }

    public static int dist(Feature2D loop) {
        return (int) Math.min(Math.abs(loop.getEnd1() - loop.getStart2()),
                Math.abs(loop.getMidPt1() - loop.getMidPt2()));
    }

}
