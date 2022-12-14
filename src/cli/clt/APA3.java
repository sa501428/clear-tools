package cli.clt;

import javastraw.feature2D.Feature2D;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.ArrayList;
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
    }
}
