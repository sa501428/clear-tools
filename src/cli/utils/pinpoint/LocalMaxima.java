package cli.utils.pinpoint;

import cli.utils.general.Pixel;

public class LocalMaxima {
    final Pixel maxCoordinate;
    final int area;

    public LocalMaxima(Pixel maxCoordinate, int area) {
        this.maxCoordinate = maxCoordinate;
        this.area = area;
    }
}
