package cli.utils.apa;

import javastraw.reader.basics.Chromosome;

public class AnchorAPAScore {

    private final Chromosome chromosome;
    private final int x1, x2;
    private final String id;
    private final float score;
    private final boolean isUpstream;

    public AnchorAPAScore(Chromosome chromosome, int x1, int x2, String id, float score, boolean isUpstream) {
        this.chromosome = chromosome;
        this.x1 = x1;
        this.x2 = x2;
        this.id = id;
        this.score = score;
        this.isUpstream = isUpstream;
    }

    public String getLineForBedgraphFile() {
        return chromosome.getName() + "\t" + x1 + "\t" + x2 + "\t" + score;
    }

    public String getLineForBEDFile() {
        return chromosome.getName() + "\t" + x1 + "\t" + x2 + "\t" +
                id + "\t" + score + "\t" + (isUpstream ? "+" : "-") +
                x1 + "\t" + x2 + "\t" + getColor();
    }

    public String getColor() {
        if (isUpstream) {
            return "0,255,0";
        } else {
            return "255,0,0";
        }
    }
}
