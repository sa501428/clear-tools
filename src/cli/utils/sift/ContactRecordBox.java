package cli.utils.sift;

import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.awt.*;
import java.util.HashMap;
import java.util.Objects;

public class ContactRecordBox {

    private final ContactRecord record;
    private final int resolution;

    public ContactRecordBox(ContactRecord record, int resolution) {
        this.record = record;
        this.resolution = resolution;
    }

    public int getResolution() {
        return resolution;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        ContactRecordBox o2 = (ContactRecordBox) o;
        return record.getBinX() == o2.record.getBinX()
                && record.getBinY() == o2.record.getBinY()
                && record.getCounts() == o2.record.getCounts()
                && resolution == o2.resolution;
    }

    @Override
    public int hashCode() {
        return Objects.hash(record.getBinX(), record.getBinY(), record.getCounts(), resolution);
    }

    public int getGenomeX1() {
        return record.getBinX() * resolution - resolution;
    }

    public int getGenomeX2() {
        return record.getBinX() * resolution + 2 * resolution;
    }

    public int getGenomeY1() {
        return record.getBinY() * resolution - resolution;
    }

    public int getGenomeY2() {
        return record.getBinY() * resolution + 2 * resolution;
    }

    public float getCounts() {
        return record.getCounts();
    }

    public Feature2D toFeature2D(Chromosome c1) {
        return new Feature2D(Feature2D.FeatureType.PEAK, c1.getName(), getGenomeX1(), getGenomeX2(),
                c1.getName(), getGenomeY1(), getGenomeY2(), Color.BLACK, new HashMap<>());
    }
}
