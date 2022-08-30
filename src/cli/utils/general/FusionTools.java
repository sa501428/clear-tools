package cli.utils.general;

import cli.clt.Cleaner;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;

public class FusionTools {

    private static final int MAX_RADIUS = 1000;

    public static void coalesceFeatures(String[] fileNames, String genomeID, String outFile, boolean useNMS) {
        Feature2DList list = combineAll(fileNames, genomeID);
        if (useNMS) {
            list.filterLists((chr, feature2DList) -> removeOverlappingPixels(feature2DList));
        } else {
            list.filterLists((chr, feature2DList) -> coalescePixelsToCentroid(feature2DList));
        }
        list.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static List<Feature2D> removeOverlappingPixels(List<Feature2D> features) {
        Map<SimpleLocation, LinkedList<Feature2D>> map = groupNearbyRecords(simpleFilter(features), 1000000);
        Set<Feature2D> coalesced = new HashSet<>();

        for (LinkedList<Feature2D> featureLL : map.values()) {

            featureLL.sort((o1, o2) -> {
                int val = Long.compare(o1.getWidth1(), o2.getWidth1());
                if (val == 0) {
                    return Long.compare(o1.getWidth2(), o2.getWidth2());
                }
                return val;
            });

            while (!featureLL.isEmpty()) {

                Feature2D pixel = featureLL.pollFirst();
                if (pixel != null) {
                    featureLL.remove(pixel);
                    Set<Feature2D> pixelList = new HashSet<>();
                    pixelList.add(pixel);

                    for (Feature2D px : featureLL) {
                        if (hasOverlap(px, pixel)) {
                            pixelList.add(px);
                        }
                    }
                    featureLL.removeAll(pixelList);
                    pixel.addStringAttribute("NumCollapsed", String.valueOf(pixelList.size()));
                    coalesced.add(pixel);
                }
            }
        }
        map.clear();

        return new ArrayList<>(coalesced);
    }

    public static List<Feature2D> coalescePixelsToCentroid(List<Feature2D> features) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        Map<SimpleLocation, LinkedList<Feature2D>> map = groupNearbyRecords(simpleFilter(features), 1000000);
        Set<Feature2D> coalesced = new HashSet<>();

        for (LinkedList<Feature2D> featureLL : map.values()) {

            long clusterRadius = Math.min(getClusterWidth(features), MAX_RADIUS);

            while (!featureLL.isEmpty()) {

                Feature2D pixel = featureLL.pollFirst();
                if (pixel != null) {
                    featureLL.remove(pixel);
                    Set<Feature2D> pixelList = new HashSet<>();
                    pixelList.add(pixel);

                    DescriptiveStatistics stats1 = getStats(pixelList, true);
                    DescriptiveStatistics stats2 = getStats(pixelList, false);

                    long median1 = (long) stats1.getPercentile(50);
                    long median2 = (long) stats2.getPercentile(50);

                    int prevSize = 0;

                    while (prevSize != pixelList.size()) {
                        prevSize = pixelList.size();
                        for (Feature2D px : featureLL) {
                            if (distance(median1 - px.getMidPt1(),
                                    median2 - px.getMidPt2()) <= clusterRadius) {
                                pixelList.add(px);
                            }
                        }
                        featureLL.removeAll(pixelList);
                        stats1 = getStats(pixelList, true);
                        stats2 = getStats(pixelList, false);
                        median1 = (long) stats1.getPercentile(50);
                        median2 = (long) stats2.getPercentile(50);
                    }

                    long r = getRadius(pixelList, median1, median2);
                    long halfWidth = getHalfWidth(pixelList);

                    long start1 = (long) (stats1.getMin() - halfWidth);
                    long start2 = (long) (stats2.getMin() - halfWidth);
                    long end1 = (long) (stats1.getMax() + halfWidth);
                    long end2 = (long) (stats2.getMax() + halfWidth);

                    Map<String, String> attributes = new HashMap<>(4);
                    attributes.put("Radius", String.valueOf(r));
                    attributes.put("Centroid1", String.valueOf(median1));
                    attributes.put("Centroid2", String.valueOf(median2));
                    attributes.put("NumCollapsed", String.valueOf(pixelList.size()));

                    coalesced.add(new Feature2D(Feature2D.FeatureType.PEAK, pixel.getChr1(), start1, end1,
                            pixel.getChr2(), start2, end2, Color.BLACK, attributes));
                }
            }
        }
        map.clear();

        return new ArrayList<>(coalesced);
    }

    private static Set<Feature2D> simpleFilter(List<Feature2D> features) {
        Set<Feature2D> featureSet = new HashSet<>();
        for (Feature2D feature : features) {
            if (Cleaner.passesMinLoopSize(feature)) {
                featureSet.add(feature);
            }
        }
        return featureSet;
    }

    private static long getHalfWidth(Set<Feature2D> pixelList) {
        long width = 0;
        for (Feature2D px2 : pixelList) {
            width = Math.max(width, px2.getWidth1());
            width = Math.max(width, px2.getWidth2());
        }
        return width / 2;
    }

    private static long getClusterWidth(List<Feature2D> pixelList) {
        long width = 0;
        for (Feature2D px2 : pixelList) {
            width = Math.max(width, px2.getWidth1());
            width = Math.max(width, px2.getWidth2());
        }
        return (long) (width * 1.5);
    }

    private static long getRadius(Set<Feature2D> pixelList, long pixelListMid1, long pixelListMid2) {
        long r = 0;
        for (Feature2D px2 : pixelList) {
            long dist = distance(pixelListMid1 - px2.getMidPt1(), pixelListMid2 - px2.getMidPt2());
            r = Math.max(r, dist);
        }
        return r;
    }

    private static DescriptiveStatistics getStats(Set<Feature2D> pixelList, boolean is1) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (Feature2D px : pixelList) {
            if (is1) {
                stats.addValue(px.getMidPt1());
            } else {
                stats.addValue(px.getMidPt2());
            }
        }
        return stats;
    }

    public static long distance(long x, long y) {
        return (long) Math.sqrt(x * x + y * y);
    }


    private static boolean hasOverlap(Feature2D px1, Feature2D px2) {
        return getWidth(px1.getStart1(), px1.getEnd1(), px2.getStart1(), px2.getEnd1()) *
                getWidth(px1.getStart2(), px1.getEnd2(), px2.getStart2(), px2.getEnd2()) > 0;
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

    public static Feature2DList combineAll(String[] fileNames, String genomeID) {
        final Feature2DList combinedLoops = new Feature2DList();

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        for (String path : fileNames) {
            Feature2DList loopList = Feature2DParser.loadFeatures(path, handler,
                    false, null, false);
            loopList.processLists(combinedLoops::addByKey);
        }

        return combinedLoops;
    }
}
