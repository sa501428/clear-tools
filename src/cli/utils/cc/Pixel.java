package cli.utils.cc;

import cli.utils.general.Utils;

public class Pixel {
    public final int row;
    public final int col;
    public final float value;
    public final float zScore;

    Pixel(int row, int col, float value, float zScore) {
        this.row = row;
        this.col = col;
        this.value = value;
        this.zScore = zScore;
    }

    public static boolean contains(Pixel px, int minR, int maxR, int minC, int maxC) {
        return Utils.inBounds(px.row, minR, maxR) && Utils.inBounds(px.col, minC, maxC);
    }
}
