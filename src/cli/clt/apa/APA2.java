package cli.clt.apa;

import cli.clt.CommandLineParser;
import cli.utils.data.BoundsInfo;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class APA2 extends APA {

    public APA2(String[] args, CommandLineParser parser) {
        super(args, parser, true);
    }

    @Override
    protected void processLoopsForRegion(MatrixZoomData zd, List<Feature2D> loops,
                                         float[][] output, AtomicInteger currNumLoops,
                                         int numTotalLoops) {
        Map<Integer, List<BoundsInfo>> loopListAsMap = convertToMap(loops);
        Set<Integer> allColumns = getColumns(loops);
        int counter = 0;

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (loopListAsMap.containsKey(cr.getBinX())) {
                    if (allColumns.contains(cr.getBinY())) {
                        populateMatrixIfApplicable(output, cr, loopListAsMap.get(cr.getBinX()));
                        if (counter++ % 1000 == 0) {
                            System.out.print(".");
                        }
                    }
                }
            }
        }

        currNumLoops.addAndGet(loops.size());

        loopListAsMap.clear();
        allColumns.clear();
        loopListAsMap = null;
        allColumns = null;
    }

    private void populateMatrixIfApplicable(float[][] matrix, ContactRecord cr, List<BoundsInfo> allBounds) {
        for (BoundsInfo bound : allBounds) {
            if (bound.contains(cr.getBinY())) {
                int relativeX = cr.getBinX() - bound.getBinXStart();
                int relativeY = cr.getBinY() - bound.getBinYStart();
                matrix[relativeX][relativeY] += cr.getCounts();
            }
        }
    }

    private Set<Integer> getColumns(List<Feature2D> loops) {
        Set<Integer> columns = new HashSet<>();
        for (Feature2D loop : loops) {
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binYEnd = binYStart + matrixWidthL;

            for (int c = binYStart; c < binYEnd; c++) {
                columns.add(c);
            }
        }
        return columns;
    }

    private Map<Integer, List<BoundsInfo>> convertToMap(List<Feature2D> loops) {
        Map<Integer, List<BoundsInfo>> loopListAsMap = new HashMap<>();
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binXEnd = binXStart + matrixWidthL;
            int binYEnd = binYStart + matrixWidthL;

            for (int r = binXStart; r < binXEnd; r++) {
                if (!loopListAsMap.containsKey(r)) {
                    loopListAsMap.put(r, new LinkedList<>());
                }
                loopListAsMap.get(r).add(new BoundsInfo(binYStart, binYEnd, binXStart));
            }
        }
        return loopListAsMap;
    }
}
