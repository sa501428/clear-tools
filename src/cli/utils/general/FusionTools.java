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
                                        boolean useNMS, boolean noAttributes, boolean useExact, boolean addIDs,
                                        String[] attributes, int roundVal) {
        Feature2DList list = combineAll(fileNames, genomeID, noAttributes, addIDs);
        if (roundVal > 1) {
            list.filterLists((chr, feature2DList) -> roundValues(feature2DList, roundVal));
        }
        list.filterLists((chr, feature2DList) -> removeOverlappingPixels(feature2DList, useNMS, useExact));
        if (attributes != null && attributes.length > 0) {
            list.filterLists((chr, feature2DList) -> filterAttributes(feature2DList, attributes));
        }
        list.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> roundValues(List<Feature2D> features, int value) {
        List<Feature2D> rounded = new ArrayList<>();
        for (Feature2D feature : features) {
            long s1 = value * (feature.getStart1() / value);
            long e1 = Math.max(value * (feature.getEnd1() / value), s1 + value);
            long s2 = value * (feature.getStart2() / value);
            long e2 = Math.max(value * (feature.getEnd2() / value), s2 + value);
            Feature2D newFeature = new Feature2D(feature.getFeatureType(),
                    feature.getChr1(), s1, e1,
                    feature.getChr2(), s2, e2,
                    feature.getColor(), feature.getAttributes());
            rounded.add(newFeature);
        }
        return rounded;
    }

    private static List<Feature2D> filterAttributes(List<Feature2D> feature2DList, String[] attributes) {
        List<Feature2D> filtered = new ArrayList<>();
        for (Feature2D feature : feature2DList) {
            Map<String, String> filteredAttributes = getSpecificAttributes(feature.getAttributes(), attributes);
            Feature2D newFeature = new Feature2D(feature.getFeatureType(), feature.getChr1(), feature.getStart1(), feature.getEnd1(),
                    feature.getChr2(), feature.getStart2(), feature.getEnd2(), feature.getColor(),
                    filteredAttributes);
            filtered.add(newFeature);
        }
        return filtered;
    }

    private static Map<String, String> getSpecificAttributes(Map<String, String> attributes, String[] attributesToKeep) {
        Map<String, String> filteredAttributes = new HashMap<>();
        for (String attribute : attributesToKeep) {
            filteredAttributes.put(attribute, attributes.getOrDefault(attribute, "NA"));
        }
        return filteredAttributes;
    }

    private static List<Feature2D> removeOverlappingPixels(List<Feature2D> features, boolean useNMS, boolean useExact) {
        if (useExact) {
            return new ArrayList<>(new HashSet<>(features));
        } else {
            Map<SimpleLocation, List<Feature2D>> map = QuickGrouping.groupNearbyRecords(features, 5000000);
            Set<Feature2D> allCoalesced = new HashSet<>();
            List<List<Feature2D>> featureLists = new ArrayList<>(map.values());
            AtomicInteger index = new AtomicInteger(0);
            ParallelizationTools.launchParallelizedCode(() -> {
                Set<Feature2D> coalesced = new HashSet<>();
                int i = index.getAndIncrement();
                while (i < featureLists.size()) {
                    List<Feature2D> featureLL = featureLists.get(i);
                    processFeatures(coalesced, featureLL, useNMS);
                    i = index.getAndIncrement();
                }
                synchronized (allCoalesced) {
                    allCoalesced.addAll(coalesced);
                }
            });
            map.clear();
            return new ArrayList<>(allCoalesced);
        }
    }

    private static void processFeatures(Set<Feature2D> coalesced, List<Feature2D> featureLL, boolean useNMS) {
        sort(featureLL);

        while (!featureLL.isEmpty()) {
            Feature2D pixel = featureLL.get(0);
            if (pixel != null) {
                featureLL.remove(pixel);
                int buffer = getBuffer(useNMS, pixel);

                List<Feature2D> pixelList = OverlapTools.getMatchesWithOverlap(pixel, featureLL, buffer);
                featureLL.removeAll(pixelList);
                //pixel.addStringAttribute("NumCollapsed", String.valueOf(pixelList.size()));

                if (useNMS) {
                    coalesced.add(pixel);
                } else {
                    pixelList.add(pixel);
                    coalesced.add(getFeatureFromBounds(pixelList));
                }
                pixelList.clear();
            }
        }
    }

    public static Feature2D getFeatureFromBounds(List<Feature2D> pixelList) {
        Feature2D pixel = pixelList.get(pixelList.size() - 1);
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
