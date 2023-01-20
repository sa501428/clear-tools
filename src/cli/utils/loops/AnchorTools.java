package cli.utils.loops;

import cli.utils.flags.Anchor;
import cli.utils.flags.LoopGenerator;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.ParallelizationTools;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorTools {

    public static Feature2DList createLoops(ChromosomeHandler handler, GenomeWide1DList<Anchor> forwardAnchors,
                                            GenomeWide1DList<Anchor> reverseAnchors,
                                            long minDist, long maxDist, int resolution) {
        System.out.println("Number of anchors: " + forwardAnchors.size() + " - " + reverseAnchors.size());

        Feature2DList output = new Feature2DList();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = index.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chromosome = chromosomes[currIndex];
                if (forwardAnchors.size() > 0 && reverseAnchors.size() > 0) {
                    List<Feature2D> newLoops = generate(chromosome, forwardAnchors, reverseAnchors, minDist, maxDist, resolution);
                    if (newLoops.size() > 0) {
                        synchronized (output) {
                            output.addByKey(Feature2DList.getKey(chromosome, chromosome), newLoops);
                        }
                    }
                }
                currIndex = index.getAndIncrement();
            }
        });
        return output;
    }

    public static List<Feature2D> generate(Chromosome chromosome, GenomeWide1DList<Anchor> forwardAnchors,
                                           GenomeWide1DList<Anchor> reverseAnchors, long minGenomeDist,
                                           long maxGenomeDist, int resolution) {

        List<Anchor> forwards = forwardAnchors.getFeatures("" + chromosome.getIndex());
        forwards.sort(Comparator.comparingLong(Anchor::getMid));

        List<Anchor> reverses = reverseAnchors.getFeatures("" + chromosome.getIndex());
        reverses.sort(Comparator.comparingLong(Anchor::getMid));

        List<Feature2D> results = new ArrayList<>();

        for (Anchor forward : forwards) {
            for (Anchor reverse : reverses) {
                if (forward.getEnd() < reverse.getStart()) {
                    Feature2D feature2D = LoopGenerator.createIntraFeature(chromosome, forward, reverse,
                            minGenomeDist, maxGenomeDist, resolution);
                    if (feature2D != null) {
                        results.add(feature2D);
                    }
                }
            }
        }
        return results;
    }

    public static Anchor getAnchor(Feature2D feature, String startPos, String endPos, int chrIndex) {
        return new Anchor(feature.getChr1(), Long.parseLong(feature.getAttribute(startPos)),
                Long.parseLong(feature.getAttribute(endPos)), chrIndex);
    }
}
