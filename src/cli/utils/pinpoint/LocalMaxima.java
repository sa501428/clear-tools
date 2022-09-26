package cli.utils.pinpoint;

import cli.utils.general.Pixel;

public class LocalMaxima {
    final Pixel maxCoordinate;
    final int area;
    final int minR, minC, maxR, maxC;

    public LocalMaxima(Pixel maxCoordinate, int area, int minR, int minC, int maxR, int maxC) {
        this.maxCoordinate = maxCoordinate;
        this.area = area;
        this.minR = minR;
        this.minC = minC;
        this.maxR = maxR;
        this.maxC = maxC;
    }
}
