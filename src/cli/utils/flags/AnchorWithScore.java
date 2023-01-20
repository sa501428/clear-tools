package cli.utils.flags;

public class AnchorWithScore extends Anchor {
    private final float score;

    public AnchorWithScore(String chromName, int start1, int end1, float score, int chrIndex) {
        super(chromName, start1, end1, chrIndex);
        this.score = score;
    }

    public float getScore() {
        return score;
    }
}
