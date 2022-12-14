package cli.clt;

import cli.utils.general.Utils;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class APA3 extends APA {
    public APA3(String[] args, CommandLineParser parser) {
        super(args, parser, true);
    }

    public static List<float[][]> initList(int size) {
        List<float[][]> list = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            list.add(new float[size][size]);
        }
        return list;
    }

    @Override
    protected void processLoopsForRegion(MatrixZoomData zd, List<Feature2D> loops,
                                         float[][] output, AtomicInteger currNumLoops,
                                         int numTotalLoops) {

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {

            }
        }

        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            Utils.addLocalBoundedRegion(output, zd, binXStart, binYStart, matrixWidthL, norm);
            if (currNumLoops.incrementAndGet() % 100 == 0) {
                System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
            }
        }
    }
}
