package cli.utils.sift;

import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.awt.*;
import java.util.HashMap;
import java.util.Objects;

public class Region {
    private static final int MIN_RESOLUTION = 100;
    private final int genomeX, genomeY, resolution;

    public Region(int binX, int binY, int resolution) {
        this.genomeX = binX * resolution;
        this.genomeY = binY * resolution;
        this.resolution = resolution;
    }

    public Region(ContactRecord cr, int resolution) {
        this(cr.getBinX(), cr.getBinY(), resolution);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        Region o2 = (Region) o;
        return genomeX == o2.genomeX && genomeY == o2.genomeY && resolution == o2.resolution;
    }

    @Override
    public int hashCode() {
        return Objects.hash(genomeX, genomeY, resolution);
    }

    public int getGenomeX() {
        return this.genomeX;
    }

    public int getGenomeY() {
        return this.genomeY;
    }

    public int getResolution() {
        return resolution;
    }

    public Feature2D toFeature2D(Chromosome c1) {
        return new Feature2D(Feature2D.FeatureType.PEAK, c1.getName(), genomeX, genomeX + resolution,
                c1.getName(), genomeY, genomeY + resolution, Color.BLACK, new HashMap<>());
    }

    public boolean containedBy(SimpleLocation location, int res) {
        for (int gx = location.getBinX() * res; gx < location.getBinX() * res + res; gx += MIN_RESOLUTION) {
            for (int gy = location.getBinY() * res; gy < location.getBinY() * res + res; gy += MIN_RESOLUTION) {
                if (gx == genomeX && gy == genomeY) {
                    return true;
                }
            }
        }
        return false;
    }
}
