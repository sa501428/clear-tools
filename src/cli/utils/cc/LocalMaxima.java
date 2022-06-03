package cli.utils.cc;

public class LocalMaxima {
    final Pixel maxCoordinate;
    final float maxVal;
    final int area;

    public LocalMaxima(Pixel maxCoordinate, float maxVal, int area) {
        this.maxCoordinate = maxCoordinate;
        this.maxVal = maxVal;
        this.area = area;
    }
}
