package cli.utils.general;

import java.util.List;

public class Pixel {
    public final int row;
    public final int col;
    public final float value;
    public final float zScore;

    public Pixel(int row, int col, float value, float zScore) {
        this.row = row;
        this.col = col;
        this.value = value;
        this.zScore = zScore;
    }

    public static boolean contains(Pixel px, int minR, int maxR, int minC, int maxC) {
        return Utils.inBounds(px.row, minR, maxR) && Utils.inBounds(px.col, minC, maxC);
    }

    public static Pixel getMax(List<Pixel> pixels) {
        Pixel max = pixels.get(0);
        for (Pixel px : pixels) {
            if (px.value > max.value) {
                max = px;
            }
        }
        return max;
    }
}
