package cli.utils.general;

import cli.utils.FeatureStats;
import cli.utils.clean.LoopTools;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class FusionTools {

    public static void coalesceFeatures(String[] fileNames, String genomeID, String outFile,
                                        boolean useNMS, boolean noAttributes, boolean useExact) {
        Feature2DList list = combineAll(fileNames, genomeID, noAttributes);
        list.filterLists((chr, feature2DList) -> removeOverlappingPixels(feature2DList, useNMS, useExact));
        list.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> removeOverlappingPixels(List<Feature2D> features, boolean useNMS, boolean useExact) {
        Map<SimpleLocation, LinkedList<Feature2D>> map = groupNearbyRecords(new HashSet<>(features), 5000000);
        Set<Feature2D> allCoalesced = new HashSet<>();

        List<LinkedList<Feature2D>> featureLists = new ArrayList<>(map.values());

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {

            Set<Feature2D> coalesced = new HashSet<>();
            int i = index.getAndIncrement();
            while (i < featureLists.size()) {
                LinkedList<Feature2D> featureLL = featureLists.get(i);
                processFeatures(coalesced, featureLL, useNMS, useExact);
                i = index.getAndIncrement();
            }

            synchronized (allCoalesced) {
                allCoalesced.addAll(coalesced);
            }

        });

        map.clear();

        return new ArrayList<>(allCoalesced);
    }

    private static void processFeatures(Set<Feature2D> coalesced, LinkedList<Feature2D> featureLL, boolean useNMS, boolean useExact) {
        sort(featureLL);

        while (!featureLL.isEmpty()) {
            Feature2D pixel = featureLL.pollFirst();
            if (pixel != null) {
                featureLL.remove(pixel);
                int buffer = getBuffer(useNMS, pixel);

                Set<Feature2D> pixelList;
                if (useExact) {
                    pixelList = getExactMatches(pixel, featureLL);
                } else {
                    pixelList = getMatchesWithOverlap(pixel, featureLL, buffer);
                }

                featureLL.removeAll(pixelList);
                pixel.addStringAttribute("NumCollapsed", String.valueOf(pixelList.size()));

                if (useNMS || useExact) {
                    coalesced.add(pixel);
                } else {
                    coalesced.add(getFeatureFromBounds(pixelList, pixel));
                }
                pixelList.clear();
            }
        }
    }

    private static Feature2D getFeatureFromBounds(Set<Feature2D> pixelList, Feature2D pixel) {
        long start1 = FeatureStats.minStart1(pixelList);
        long start2 = FeatureStats.minStart2(pixelList);
        long end1 = FeatureStats.maxEnd1(pixelList);
        long end2 = FeatureStats.maxEnd2(pixelList);
        return new Feature2D(Feature2D.FeatureType.PEAK, pixel.getChr1(), start1, end1,
                pixel.getChr2(), start2, end2, Color.BLACK, pixel.getAttributes());
    }

    private static Set<Feature2D> getMatchesWithOverlap(Feature2D pixel, LinkedList<Feature2D> featureLL, int buffer) {
        Set<Feature2D> pixelList = new HashSet<>();
        pixelList.add(pixel);
        for (Feature2D px : featureLL) {
            if (hasOverlap(px, pixel, buffer)) {
                pixelList.add(px);
            }
        }
        return pixelList;
    }

    private static Set<Feature2D> getExactMatches(Feature2D pixel, LinkedList<Feature2D> featureLL) {
        Set<Feature2D> pixelList = new HashSet<>();
        pixelList.add(pixel);
        for (Feature2D px : featureLL) {
            if (isExact(px, pixel)) {
                pixelList.add(px);
            }
        }
        return pixelList;
    }

    private static int getBuffer(boolean useNMS, Feature2D pixel) {
        if (useNMS) {
            return 0;
        } else {
            return (int) Math.max(pixel.getWidth1(), pixel.getWidth2());
        }
    }

    private static void sort(LinkedList<Feature2D> featureLL) {
        Collections.sort(featureLL, (o1, o2) -> {
            int val = Long.compare(o1.getWidth1(), o2.getWidth1());
            if (val == 0) {
                return Long.compare(o1.getWidth2(), o2.getWidth2());
            }
            return val;
        });
    }

    public static long distance(long x, long y) {
        return (long) Math.sqrt(x * x + y * y);
    }

    private static boolean isExact(Feature2D px, Feature2D px2) {
        return px.getStart1() == px2.getStart1()
                && px.getStart2() == px2.getStart2()
                && px.getEnd1() == px2.getEnd1()
                && px.getEnd2() == px2.getEnd2();
    }

    private static boolean hasOverlap(Feature2D px1, Feature2D original, int buffer) {
        return getWidth(px1.getStart1(), px1.getEnd1(), original.getStart1() - buffer, original.getEnd1() + buffer) *
                getWidth(px1.getStart2(), px1.getEnd2(), original.getStart2() - buffer, original.getEnd2() + buffer) > 0;
    }

    private static int getWidth(long p1, long p2, long g1, long g2) {
        long start = Math.max(p1, g1);
        long end = Math.min(p2, g2);
        return (int) Math.max(0, end - start);
    }

    public static Map<SimpleLocation, LinkedList<Feature2D>> groupNearbyRecords(Set<Feature2D> initialPoints, int scalar) {
        Map<SimpleLocation, LinkedList<Feature2D>> locationMap = new HashMap<>();
        for (Feature2D feature : initialPoints) {
            SimpleLocation region = new SimpleLocation((int) (feature.getMidPt1() / scalar), (int) (feature.getMidPt2() / scalar));
            if (locationMap.containsKey(region)) {
                locationMap.get(region).add(feature);
            } else {
                LinkedList<Feature2D> values = new LinkedList<>();
                values.add(feature);
                locationMap.put(region, values);
            }
        }
        return locationMap;
    }

    public static Feature2DList combineAll(String[] fileNames, String genomeID, boolean noAttributes) {
        final Feature2DList combinedLoops = new Feature2DList();

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        for (String path : fileNames) {
            Feature2DList loopList = LoopTools.loadFilteredBedpe(path, handler, !noAttributes);
            loopList.processLists(combinedLoops::addByKey);
        }

        return combinedLoops;
    }
}
