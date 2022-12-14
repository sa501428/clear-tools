package cli.clt;

import java.util.ArrayList;
import java.util.List;

public class APA3 {

    public static List<float[][]> initList(int numLoops, int matrixSize) {
        List<float[][]> list = new ArrayList<>(numLoops);
        for (int i = 0; i < numLoops; i++) {
            list.add(new float[matrixSize][matrixSize]);
        }
        return list;
    }
}
