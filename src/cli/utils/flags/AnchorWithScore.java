package cli.utils.flags;

public class AnchorWithScore extends Anchor {
    private final float score;

    public AnchorWithScore(String chromName, int start1, int end1, float score) {
        super(chromName, start1, end1);
        this.score = score;
    }

    public float getScore() {
        return score;
    }
}
