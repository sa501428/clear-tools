package cli.utils.flags;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;

import java.awt.*;
import java.util.List;
import java.util.*;

public class LoopGenerator {
    public static List<Feature2D> generate(GenomeWide1DList<Anchor> anchors, Chromosome chrom1, Chromosome chrom2,
                                           long minGenomeDist, long maxGenomeDist, int res) {
        List<Feature2D> results = new ArrayList<>();
        if (chrom1.getIndex() == chrom2.getIndex()) {
            List<Anchor> aList = anchors.getFeatures("" + chrom1.getIndex());
            aList = new ArrayList<>(aList);
            aList.sort(Comparator.comparingLong(Anchor::getMid));
            for (int i = 0; i < aList.size(); i++) {
                for (int j = i + 1; j < aList.size(); j++) {
                    Feature2D feature2D = createIntraFeature(chrom1, aList.get(i), aList.get(j),
                            minGenomeDist, maxGenomeDist, res);
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

    public static Feature2D createIntraFeature(Chromosome chrom1, Anchor a1, Anchor a2,
                                               long minGenomeDist, long maxGenomeDist, int resolution) {
        long dist = a2.getMid() - a1.getMid();
        if (dist < 0) {
            System.err.println("Weird error!!");
            System.exit(10);
        }
        if (dist > minGenomeDist && dist <= maxGenomeDist) {
            Map<String, String> attributes = new HashMap<>(4);
            attributes.put("motif_start_1", "" + a1.getStart());
            attributes.put("motif_end_1", "" + a1.getEnd());
            attributes.put("motif_start_2", "" + a2.getStart());
            attributes.put("motif_end_2", "" + a2.getEnd());
            if (resolution > 1) {
                long x1 = round(a1.getStart(), resolution);
                long y1 = round(a2.getStart(), resolution);
                long x2 = Math.max(round(a1.getEnd(), resolution) + 1, x1 + resolution);
                long y2 = Math.max(round(a2.getEnd(), resolution) + 1, y1 + resolution);
                return new Feature2D(Feature2D.FeatureType.PEAK,
                        chrom1.getName(), x1, x2, chrom1.getName(), y1, y2, Color.BLACK, attributes);
            } else {
                return new Feature2D(Feature2D.FeatureType.PEAK, chrom1.getName(), a1.getStart(), a1.getEnd(),
                        chrom1.getName(), a2.getStart(), a2.getEnd(), Color.BLACK, attributes);
            }
        }
        return null;
    }

    private static long round(long number, int resolution) {
        return (number / resolution) * resolution;
    }

    private static Feature2D createFeature(Chromosome chrom1, Anchor a1, Chromosome chrom2, Anchor a2) {
        return new Feature2D(Feature2D.FeatureType.PEAK, chrom1.getName(), a1.getStart(), a1.getEnd(),
                chrom2.getName(), a2.getStart(), a2.getEnd(), Color.BLACK, new HashMap<>(0));
    }
}
