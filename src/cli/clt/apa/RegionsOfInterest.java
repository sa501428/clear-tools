package cli.clt.apa;

import cli.utils.data.Bounds2DInfo;
import cli.utils.data.CompressedGridMapForBounds2D;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Class to manage regions of interest with respect to a grid map and contact records.
 */
public class RegionsOfInterest {
    private final int maxNumForIndices;
    private final int resolution;
    private final int window;
    private final int matrixWidthL;
    private final List<CompressedGridMapForBounds2D> loopListsAsMaps;
    private final Set<Integer> allRowIndices = new HashSet<>();
    private final Set<Integer> allColIndices = new HashSet<>();

    public RegionsOfInterest(int resolution, int window, int matrixWidthL,
                             List<List<Feature2D>> allLoops) {
        maxNumForIndices = allLoops.size();
        this.resolution = resolution;
        this.window = window;
        this.matrixWidthL = matrixWidthL;
        loopListsAsMaps = convertToMaps(allLoops);
    }

    private List<CompressedGridMapForBounds2D> convertToMaps(List<List<Feature2D>> allLoops) {
        List<CompressedGridMapForBounds2D> loopListsAsMaps = new ArrayList<>(allLoops.size());
        for (List<Feature2D> loops : allLoops) {
            loopListsAsMaps.add(new CompressedGridMapForBounds2D(loops, resolution, window, matrixWidthL,
                    allRowIndices, allColIndices));
        }
        return loopListsAsMaps;
    }

    public boolean containsRecord(ContactRecord cr, int i) {
        validateIndex(i);
        return loopListsAsMaps.get(i).contains(cr.getBinX(), cr.getBinY());
    }

    public List<Bounds2DInfo> getBoundsInfo(ContactRecord cr, int i) {
        validateIndex(i);
        return loopListsAsMaps.get(i).get(cr.getBinX(), cr.getBinY());
    }

    public boolean probablyContainsRecord(ContactRecord cr) {
        return allRowIndices.contains(cr.getBinX()) && allColIndices.contains(cr.getBinY());
    }

    public void clear() {
        loopListsAsMaps.clear();
        allRowIndices.clear();
        allColIndices.clear();
    }

    private void validateIndex(int i) {
        if (i < 0 || i >= maxNumForIndices) {
            throw new IndexOutOfBoundsException("Index " + i + " is out of bounds. Valid indices are 0 to " + (loopListsAsMaps.size() - 1));
        }
    }
}
