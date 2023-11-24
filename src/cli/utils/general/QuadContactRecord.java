package cli.utils.general;

import javastraw.reader.block.ContactRecord;

public class QuadContactRecord {
    private final int binX;
    private final int binY;
    private final float counts;
    private final float oe;
    private final float perc;
    private final float zscore;


    public QuadContactRecord(ContactRecord cr, float oe, float perc, float zscore) {
        this.binX = cr.getBinX();
        this.binY = cr.getBinY();
        this.counts = cr.getCounts();
        this.oe = oe;
        this.perc = perc;
        this.zscore = zscore;
    }

    public QuadContactRecord(Integer i, Integer j, float maxCounts, float maxOE, float maxPerc, float maxZscore) {
        this.binX = i;
        this.binY = j;
        this.counts = maxCounts;
        this.oe = maxOE;
        this.perc = maxPerc;
        this.zscore = maxZscore;
    }

    public int getBinX() {
        return binX;
    }

    public int getBinY() {
        return binY;
    }

    public float getCounts() {
        return counts;
    }

    public float getOE() {
        return oe;
    }

    public float getPerc() {
        return perc;
    }

    public float getZscore() {
        return zscore;
    }
}
