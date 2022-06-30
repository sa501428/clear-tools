package cli.utils;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.Collection;

public class StandardDevUtils {
    public static float[][] StdDeviationFinder(Collection<float[][]> mtxList, int window) {
        /*
        Computes a 2D float array that stores the standard deviation of each respective pixel across all Hi-C maps in mtxList

        Inputs:
            mtxList -- a List containing 2D float arrays that represent Hi-C maps

        Outputs a 2D float array the stores the standard deviation of each respective pixel across all Hi-C maps in mtxList
        */
        float[][] output = new float[window][window];
        float[][] range = new float[window][window];
        for (int row = 0; row < window; row++) {
            for (int col = 0; col < window; col++) {
                DescriptiveStatistics stats = new DescriptiveStatistics();
                for (float[][] mtx : mtxList) {
                    stats.addValue(mtx[row][col]);
                }
                output[row][col] = (float) stats.getStandardDeviation();
                range[row][col] = (float) (stats.getMax() - stats.getMin());
            }
        }
        return output;
    }
}
