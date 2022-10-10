package cli.utils.general;

import cli.utils.FeatureStats;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
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
                                        boolean useNMS, boolean noAttributes, boolean useExact, boolean addIDs) {
        Feature2DList list = combineAll(fileNames, genomeID, noAttributes, addIDs);
        list.filterLists((chr, feature2DList) -> removeOverlappingPixels(feature2DList, useNMS, useExact));
        list.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> removeOverlappingPixels(List<Feature2D> features, boolean useNMS, boolean useExact) {
        Map<SimpleLocation, List<Feature2D>> map = QuickGrouping.groupNearbyRecords(features, 5000000);
        Set<Feature2D> allCoalesced = new HashSet<>();

        List<List<Feature2D>> featureLists = new ArrayList<>(map.values());

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {

            Set<Feature2D> coalesced = new HashSet<>();
            int i = index.getAndIncrement();
            while (i < featureLists.size()) {
                List<Feature2D> featureLL = featureLists.get(i);
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

    private static void processFeatures(Set<Feature2D> coalesced, List<Feature2D> featureLL, boolean useNMS, boolean useExact) {
        sort(featureLL);

        while (!featureLL.isEmpty()) {
            Feature2D pixel = featureLL.get(0);
            if (pixel != null) {
                featureLL.remove(pixel);
                int buffer = getBuffer(useNMS, pixel);

                Set<Feature2D> pixelList;
                if (useExact) {
                    pixelList = OverlapTools.getExactMatches(pixel, featureLL);
                } else {
                    pixelList = OverlapTools.getMatchesWithOverlap(pixel, featureLL, buffer);
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

    private static int getBuffer(boolean useNMS, Feature2D pixel) {
        if (useNMS) {
            return 0;
        } else {
            return (int) Math.max(pixel.getWidth1(), pixel.getWidth2());
        }
    }

    private static void sort(List<Feature2D> featureLL) {
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

    public static Feature2DList combineAll(String[] fileNames, String genomeID, boolean noAttributes,
                                           boolean addIDs) {
        final Feature2DList combinedLoops = new Feature2DList();

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        int counter = 0;
        for (String path : fileNames) {
            Feature2DList loopList = Feature2DParser.loadFeatures(path, handler, !noAttributes, null, false);
            if (addIDs) {
                addIDToAllLoops(loopList, counter, "LIST_UID");
                counter++;
            }
            loopList.processLists(combinedLoops::addByKey);
        }

        return combinedLoops;
    }

    private static void addIDToAllLoops(Feature2DList loopList, int id, String attribute) {
        loopList.processLists((s, list) -> {
            for (Feature2D feature : list) {
                feature.addIntAttribute(attribute, id);
            }
        });
    }
}
