package cli.utils.sift;

import cli.utils.WelfordStats;

public class Z4Scores {
    private final ZScores zScoresRaw, zScoresVC, zScoresVCSqrt, zScoresScale;

    public Z4Scores(WelfordStats statsRaw, WelfordStats statsVC, WelfordStats statsVCSqrt, WelfordStats statsScale) {
        zScoresRaw = statsRaw.getZscores();
        zScoresVC = statsVC.getZscores();
        zScoresVCSqrt = statsVCSqrt.getZscores();
        zScoresScale = statsScale.getZscores();
    }

    public boolean passesAllZscores(int index, int cutoff, double raw, double vc, double vcsqrt, double scale) {
        return zScoresRaw.getZscore(index, raw) > cutoff &&
                zScoresVC.getZscore(index, vc) > cutoff &&
                zScoresVCSqrt.getZscore(index, vcsqrt) > cutoff &&
                zScoresScale.getZscore(index, scale) > cutoff;
    }
}
