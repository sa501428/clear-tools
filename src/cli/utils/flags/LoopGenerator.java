package cli.utils.flags;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;

import java.awt.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class LoopGenerator {
    public static List<Feature2D> generate(GenomeWide1DList<Anchor> anchors, Chromosome chrom1, Chromosome chrom2,
                                           long minGenomeDist, long maxGenomeDist) {
        List<Feature2D> results = new ArrayList<>();
        if (chrom1.getIndex() == chrom2.getIndex()) {
            List<Anchor> aList = anchors.getFeatures("" + chrom1.getIndex());
            aList = new ArrayList<>(aList);
            aList.sort(Comparator.comparingInt(Anchor::getMid));
            for (int i = 0; i < aList.size(); i++) {
                for (int j = i + 1; j < aList.size(); j++) {
                    Feature2D feature2D = createIntraFeature(chrom1, aList.get(i), aList.get(j),
                            minGenomeDist, maxGenomeDist);
                    if (feature2D != null) {
                        results.add(feature2D);
                    }
                }
            }
        } else {
            for (Anchor a1 : anchors.getFeatures("" + chrom1.getIndex())) {
                for (Anchor a2 : anchors.getFeatures("" + chrom2.getIndex())) {
                    results.add(createFeature(chrom1, a1, chrom2, a2));
                }
            }
        }
        return results;
    }

    private static Feature2D createIntraFeature(Chromosome chrom1, Anchor a1, Anchor a2,
                                                long minGenomeDist, long maxGenomeDist) {
        int dist = a2.getMid() - a1.getMid();
        if (dist < 0) {
            System.err.println("Weird error!!");
            System.exit(10);
        }
        if (dist > minGenomeDist && dist <= maxGenomeDist) {
            return new Feature2D(Feature2D.FeatureType.PEAK, chrom1.getName(), a1.getStart(), a1.getEnd(),
                    chrom1.getName(), a2.getStart(), a2.getEnd(), Color.BLACK, null);
        }
        return null;
    }

    private static Feature2D createFeature(Chromosome chrom1, Anchor a1, Chromosome chrom2, Anchor a2) {
        return new Feature2D(Feature2D.FeatureType.PEAK, chrom1.getName(), a1.getStart(), a1.getEnd(),
                chrom2.getName(), a2.getStart(), a2.getEnd(), Color.BLACK, null);
    }
}
