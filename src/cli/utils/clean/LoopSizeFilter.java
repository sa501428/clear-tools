package cli.utils.clean;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DFilter;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;

import java.util.ArrayList;
import java.util.List;

public class LoopSizeFilter {

    private static final int MIN_LOOP_SIZE = 30000;

    public static Feature2DList loadFilteredBedpe(String bedpeFile, ChromosomeHandler handler,
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
        return floorDist(loop) > MIN_LOOP_SIZE;
    }

    public static long floorDist(Feature2D loop) {
        return Math.min(Math.abs(loop.getEnd1() - loop.getStart2()),
                Math.abs(loop.getMidPt1() - loop.getMidPt2()));
    }

}
