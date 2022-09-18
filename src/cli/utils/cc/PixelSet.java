package cli.utils.cc;

import java.util.ArrayList;
import java.util.List;

public class PixelSet {
    int x, y;
    List<Pixel> pixels = new ArrayList<>();

    public PixelSet(Pixel pixel) {
        add(pixel);
    }

    public int size() {
        return pixels.size();
    }

    public void add(Pixel px) {
        pixels.add(px);
        updateCentroid();
    }

    private void updateCentroid() {
        int xs = 0, ys = 0;
        for (Pixel p : pixels) {
            xs += p.x;
            ys += p.y;
        }
        x = xs / pixels.size();
        y = ys / pixels.size();
    }

    public int distanceTo(Pixel px) {
        return Math.abs(px.x - x) + Math.abs(px.y - y);
    }

    public Pixel toPixel() {
        return new Pixel(x, y);
    }
}
