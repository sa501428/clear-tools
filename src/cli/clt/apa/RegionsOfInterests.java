package cli.clt.apa;

import cli.utils.data.BoundsInfo;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.util.*;

public class RegionsOfInterests {
    private final int resolution;
    private final int window;
    private final int matrixWidthL;
    private final List<Map<Integer, List<BoundsInfo>>> loopListsAsMaps;
    private final Set<Integer> allRowIndices = new HashSet<>();
    private final Set<Integer> allColIndices = new HashSet<>();

    public RegionsOfInterests(int resolution, int window, int matrixWidthL, List<List<Feature2D>> allLoops) {
        this.resolution = resolution;
        this.window = window;
        this.matrixWidthL = matrixWidthL;
        loopListsAsMaps = convertToMaps(allLoops);
    }

    private List<Map<Integer, List<BoundsInfo>>> convertToMaps(List<List<Feature2D>> allLoops) {
        List<Map<Integer, List<BoundsInfo>>> loopListsAsMaps = new ArrayList<>(allLoops.size());
        for (List<Feature2D> loops : allLoops) {
            loopListsAsMaps.add(convertToMap(loops));
        }
        return loopListsAsMaps;
    }

    private Map<Integer, List<BoundsInfo>> convertToMap(List<Feature2D> loops) {
        Map<Integer, List<BoundsInfo>> loopListAsMap = new HashMap<>();
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binXEnd = binXStart + matrixWidthL;
            int binYEnd = binYStart + matrixWidthL;

            for (int r = binXStart; r < binXEnd; r++) {
                allRowIndices.add(r);
                if (!loopListAsMap.containsKey(r)) {
                    loopListAsMap.put(r, new LinkedList<>());
                }
                loopListAsMap.get(r).add(new BoundsInfo(binYStart, binYEnd, binXStart));
            }

            for (int c = binYStart; c < binYEnd; c++) {
                allColIndices.add(c);
            }
        }
        return loopListAsMap;
    }

    public boolean containsRecord(ContactRecord cr, int i) {
        return loopListsAsMaps.get(i).containsKey(cr.getBinX());
    }

    public List<BoundsInfo> getBoundsInfo(ContactRecord cr, int i) {
        return loopListsAsMaps.get(i).get(cr.getBinX());
    }

    public boolean probablyContainsRecord(ContactRecord cr) {
        return allRowIndices.contains(cr.getBinX()) && allColIndices.contains(cr.getBinY());
    }


    public void clear() {
        loopListsAsMaps.clear();
        allRowIndices.clear();
        allColIndices.clear();
    }
}
