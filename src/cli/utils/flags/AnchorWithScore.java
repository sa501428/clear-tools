package cli.utils.flags;

import javastraw.feature1D.Feature1D;

import java.util.Objects;

public class AnchorWithScore extends Anchor {
    protected final String name;
    protected final float score;

    public AnchorWithScore(String chromName, int start1, int end1, float score, int chrIndex, String name) {
        super(chromName, start1, end1, chrIndex);
        this.score = score;
        this.name = name;
    }

    public float getScore() {
        return score;
    }

    public String getName() {
        return name;
    }

    @Override
    public Feature1D deepClone() {
        return new AnchorWithScore(chromName, (int) start, (int) end, score, chrIndex, name);
    }

    @Override
    public String toString() {
        return chromName + "\t" + start + "\t" + end
                + "\t" + name + "\t" + score;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof AnchorWithScore) {
            AnchorWithScore a = (AnchorWithScore) obj;
            return a.chromName.equals(chromName) && a.start == start && a.end == end
                    && a.name.equals(name) && a.score == score;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromName, start, end, score, name);
    }
}
