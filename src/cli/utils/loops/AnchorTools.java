package cli.utils.loops;

import cli.utils.flags.Anchor;
import cli.utils.flags.LoopGenerator;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.ParallelizationTools;

import java.util.*;
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

    public static Map<Integer, List<Anchor>> getAnchorMap(List<Anchor> anchors, int compression) {
        Map<Integer, List<Anchor>> anchorMap = new HashMap<>();
        for (Anchor anchor : anchors) {
            int key = (int) (anchor.getMid() / compression);
            if (!anchorMap.containsKey(key)) {
                anchorMap.put(key, new ArrayList<>());
            }
            anchorMap.get(key).add(anchor);
        }
        return anchorMap;
    }

    public static List<Anchor> getClosestAnchors(Map<Integer, List<Anchor>> forwardsMap, long start, long end, int compression) {
        int key = (int) (((start + end) / 2) / compression);
        List<Anchor> anchors = new ArrayList<>();
        if (forwardsMap.containsKey(key)) {
            anchors.addAll(forwardsMap.get(key));
        }
        if (forwardsMap.containsKey(key - 1)) {
            anchors.addAll(forwardsMap.get(key - 1));
        }
        if (forwardsMap.containsKey(key + 1)) {
            anchors.addAll(forwardsMap.get(key + 1));
        }
        return anchorsContainedWithinRegion(start, end, anchors);
    }

    private static List<Anchor> anchorsContainedWithinRegion(long start, long end, List<Anchor> anchors) {
        List<Anchor> containedAnchors = new ArrayList<>();
        for (Anchor anchor : anchors) {
            if (anchor.getMid() >= start && anchor.getEnd() < end) {
                containedAnchors.add(anchor);
            }
        }
        return containedAnchors;
    }
}
