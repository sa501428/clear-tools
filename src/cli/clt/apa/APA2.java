package cli.clt.apa;

import cli.clt.CommandLineParser;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class APA2 extends APA {

    private final static int C0 = 0, C1 = 1, R0 = 2;

    public APA2(String[] args, CommandLineParser parser) {
        super(args, parser, true);
    }

    @Override
    protected void processLoopsForRegion(MatrixZoomData zd, List<Feature2D> loops,
                                         float[][] output, AtomicInteger currNumLoops,
                                         int numTotalLoops) {
        Map<Integer, List<int[]>> loopListAsMap = convertToMap(loops);
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

    private void populateMatrixIfApplicable(float[][] matrix, ContactRecord cr, List<int[]> allBounds) {
        for (int[] bounds : allBounds) {
            if (cr.getBinY() >= bounds[C0] && cr.getBinY() < bounds[C1]) {
                int relativeX = cr.getBinX() - bounds[R0];
                int relativeY = cr.getBinY() - bounds[C0];
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

    private Map<Integer, List<int[]>> convertToMap(List<Feature2D> loops) {
        Map<Integer, List<int[]>> loopListAsMap = new HashMap<>();
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binXEnd = binXStart + matrixWidthL;
            int binYEnd = binYStart + matrixWidthL;

            for (int r = binXStart; r < binXEnd; r++) {
                if (!loopListAsMap.containsKey(r)) {
                    loopListAsMap.put(r, new LinkedList<>());
                }
                // C0 = 0, C1 = 1, R0 = 2
                loopListAsMap.get(r).add(new int[]{binYStart, binYEnd, binXStart});
            }
        }
        return loopListAsMap;
    }
}
