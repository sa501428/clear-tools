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
                                            long minDist, long maxDist, int resolution,
                                            boolean makeIntra, boolean makeInter) {
        System.out.println("Number of anchors: " + forwardAnchors.size() + " - " + reverseAnchors.size());

        Feature2DList output = new Feature2DList();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = index.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chromosome = chromosomes[currIndex];
                if (makeIntra && forwardAnchors.size() > 0 && reverseAnchors.size() > 0) {
                    List<Feature2D> newLoops = generate(chromosome, forwardAnchors, reverseAnchors, minDist, maxDist, resolution);
                    if (newLoops.size() > 0) {
                        synchronized (output) {
                            output.addByKey(Feature2DList.getKey(chromosome, chromosome), newLoops);
                        }
                    }
                }
                if (makeInter) {
                    for (int j = currIndex + 1; j < chromosomes.length; j++) {
                        Chromosome chromosome2 = chromosomes[j];
                        List<Feature2D> newLoops = generateForInter(chromosome, chromosome2,
                                forwardAnchors, reverseAnchors);
                        if (newLoops.size() > 0) {
                            synchronized (output) {
                                output.addByKey(Feature2DList.getKey(chromosome, chromosome2), newLoops);
                            }
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

        List<Anchor> forwards = getSortedAnchors(forwardAnchors, "" + chromosome.getIndex());
        List<Anchor> reverses = getSortedAnchors(reverseAnchors, "" + chromosome.getIndex());

        List<Feature2D> results = new LinkedList<>();

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

    public static List<Feature2D> generateForInter(Chromosome chromosome1, Chromosome chromosome2,
                                                   GenomeWide1DList<Anchor> forwardAnchors,
                                                   GenomeWide1DList<Anchor> reverseAnchors) {

        List<Anchor> forwards = getSortedAnchors(forwardAnchors, "" + chromosome1.getIndex());
        List<Anchor> reverses = getSortedAnchors(reverseAnchors, "" + chromosome2.getIndex());

        List<Feature2D> results = new LinkedList<>();

        for (Anchor forward : forwards) {
            for (Anchor reverse : reverses) {
                results.add(LoopGenerator.createFeature(chromosome1, forward, chromosome2, reverse));
            }
        }
        return results;
    }

    private static List<Anchor> getSortedAnchors(GenomeWide1DList<Anchor> allAnchors, String key) {
        List<Anchor> anchors = allAnchors.getFeatures(key);
        anchors.sort(Comparator.comparingLong(Anchor::getMid));
        return anchors;
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

    public static Anchor getClosestAnchor(Map<Integer, List<Anchor>> anchorMap, long start, long end, int compression) {
        int key = (int) (((start + end) / 2) / compression);
        List<Anchor> anchors = new ArrayList<>();
        if (anchorMap.containsKey(key)) {
            anchors.addAll(anchorMap.get(key));
        }
        if (anchorMap.containsKey(key - 1)) {
            anchors.addAll(anchorMap.get(key - 1));
        }
        if (anchorMap.containsKey(key + 1)) {
            anchors.addAll(anchorMap.get(key + 1));
        }
        return closestAnchorToMid(start, end, anchors);
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

    private static Anchor closestAnchorToMid(long start, long end, List<Anchor> anchors) {
        long mid = (start + end) / 2;
        long minDist = Long.MAX_VALUE;
        Anchor closestAnchor = null;
        for (Anchor anchor : anchors) {
            long dist = Math.abs(mid - anchor.getMid());
            if (dist < minDist) {
                minDist = dist;
                closestAnchor = anchor;
            }
        }
        return closestAnchor;
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
