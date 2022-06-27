package cli.utils;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.List;

public class StandardDevUtils {
    public static float[][] StdDeviationFinder(List<float[][]> mtxList) {
        /*
        Computes a 2D float array that stores the standard deviation of each respective pixel across all Hi-C maps in mtxList

        Inputs:
            mtxList -- a List containing 2D float arrays that represent Hi-C maps

        Outputs a 2D float array the stores the standard deviation of each respective pixel across all Hi-C maps in mtxList
        */
        int numRows = mtxList.get(0).length;
        int numCols = mtxList.get(0)[0].length;
        float[][] output = new float[numRows][numCols];
        float[][] range = new float[numRows][numCols];
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
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
