package cli.utils.cc;

import javastraw.feature2D.Feature2D;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class Location2D {

    private final int area;
    private final float maxVal;
    private final String chr1, chr2;
    private final long start1, start2, end1, end2;

    public Location2D(String chr1, String chr2, long binXStart, long binYStart, LocalMaxima max, int resolution) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        start1 = resolution * (binXStart + max.maxCoordinate.x);
        start2 = resolution * (binYStart + max.maxCoordinate.y);
        end1 = start1 + resolution;
        end2 = start2 + resolution;
        this.area = max.area;
        this.maxVal = max.maxVal;
    }

    public Feature2D toFeature2D() {
        Map<String, String> attributes = new HashMap<>();
        attributes.put("pinpoint_area", "" + area);
        attributes.put("pinpoint_enrichment", "" + maxVal);
        return new Feature2D(Feature2D.FeatureType.PEAK, chr1, start1, end1,
                chr2, start2, end2, Color.BLACK, attributes);
    }

    @Override
    public int hashCode() {
        return Objects.hash(this.chr1, this.start1, this.end1, this.chr2, this.start2, this.end2);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        } else if (this.getClass() != obj.getClass()) {
            return false;
        } else if (this == obj) {
            return true;
        } else {
            Location2D other = (Location2D) obj;
            return this.chr1.equals(other.chr1) && this.chr2.equals(other.chr2)
                    && this.start1 == other.start1 && this.start2 == other.start2
                    && this.end1 == other.end1 && this.end2 == other.end2;
        }
    }
}
