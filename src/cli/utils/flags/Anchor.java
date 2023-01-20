package cli.utils.flags;

import javastraw.feature1D.Feature1D;

import java.util.Objects;

public class Anchor extends Feature1D {

    private final String chromName;
    private final long start, end, mid;

    public Anchor(String chromName, long start1, long end1) {
        if (chromName.startsWith("chr")) {
            this.chromName = chromName;
        } else {
            this.chromName = "chr" + chromName;
        }
        this.start = start1;
        this.end = end1;
        this.mid = (start1 + end1) / 2;
    }

    @Override
    public String getKey() {
        return chromName;
    }

    @Override
    public Feature1D deepClone() {
        return new Anchor(chromName, start, end);
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    public long getMid() {
        return mid;
    }

    @Override
    public String toString() {
        return chromName + "\t" + start + "\t" + end;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Anchor) {
            Anchor a = (Anchor) obj;
            return a.chromName.equals(chromName) && a.start == start && a.end == end;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromName, start, end);
    }
}
