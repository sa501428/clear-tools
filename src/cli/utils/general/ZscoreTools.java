package cli.utils.general;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;

public class ZscoreTools {

    public static double getLLZscore(float[][] regionMatrix, int midX, int midY, int window) {
        int startR = Math.max(midX + 1, 0);
        int endR = Math.min(midX + window + 1, regionMatrix.length);
        int startC = Math.max(midY - window, 0);
        int endC = Math.min(midY - 1, regionMatrix[0].length);
        return getZscoreForRegion(regionMatrix, midX, midY, startR, endR, startC, endC);
    }

    public static double getLocalZscore(float[][] regionMatrix, int midX, int midY, int window) {
        int startR = Math.max(midX - window, 0);
        int endR = Math.min(midX + window + 1, regionMatrix.length);
        int startC = Math.max(midY - window, 0);
        int endC = Math.min(midY + window + 1, regionMatrix[0].length);
        return getZscoreForRegion(regionMatrix, midX, midY, startR, endR, startC, endC);
    }

    public static double getZscoreForRegion(float[][] regionMatrix, int midX, int midY, int startR, int endR, int startC, int endC) {
        Welford welford = new Welford();
        for (int i = startR; i < endR; i++) {
            for (int j = startC; j < endC; j++) {
                if (i != midX && j != midY) {
                    welford.addValue(regionMatrix[i][j]);
                }
            }
        }
        Zscore zscore = welford.getZscore();
        return zscore.getZscore(regionMatrix[midX][midY]);
    }
}
