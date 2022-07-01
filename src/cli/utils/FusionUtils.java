package cli.utils;

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

public class FusionUtils {

    public static void coalesceFeaturesToCentroid(String[] fileNames, String genomeID, String outFile) {
        Feature2DList list = combineAll(fileNames, genomeID);
        list.filterLists((chr, feature2DList) -> coalescePixelsToCentroid(feature2DList));
        list.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
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

    public static List<Feature2D> coalescePixelsToCentroid(List<Feature2D> features) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        LinkedList<Feature2D> featureLL = new LinkedList<>(new HashSet<>(features));
        List<Feature2D> coalesced = new ArrayList<>();

        long clusterRadius = getClusterWidth(features);

        while (!featureLL.isEmpty()) {

            Feature2D pixel = featureLL.pollFirst();
            if (pixel != null) {
                featureLL.remove(pixel);
                Set<Feature2D> pixelList = new HashSet<>();
                pixelList.add(pixel);

                long pixelListMid1 = pixel.getMidPt1();
                long pixelListMid2 = pixel.getMidPt2();

                int prevSize = 0;

                while (prevSize != pixelList.size()) {
                    prevSize = pixelList.size();
                    for (Feature2D px : featureLL) {
                        if (distance(pixelListMid1 - px.getMidPt1(),
                                pixelListMid2 - px.getMidPt2()) <= clusterRadius) {
                            pixelList.add(px);
                        }
                    }
                    pixelListMid1 = median(pixelList, true);
                    pixelListMid2 = median(pixelList, false);
                    featureLL.removeAll(pixelList);
                }

                long r = getRadius(pixelList, pixelListMid1, pixelListMid2);
                long halfWidth = getHalfWidth(pixelList);

                long start1 = pixelListMid1 - halfWidth;
                long start2 = pixelListMid2 - halfWidth;
                long end1 = pixelListMid1 + halfWidth;
                long end2 = pixelListMid2 + halfWidth;

                Map<String, String> attributes = new HashMap<>(4);
                attributes.put("Radius", String.valueOf(r));
                attributes.put("Centroid1", String.valueOf(pixelListMid1));
                attributes.put("Centroid2", String.valueOf(pixelListMid2));
                attributes.put("NumCollapsed", String.valueOf(pixelList.size()));

                coalesced.add(new Feature2D(Feature2D.FeatureType.PEAK, pixel.getChr1(), start1, end1,
                        pixel.getChr2(), start2, end2, Color.BLACK, attributes));
            }
        }

        return coalesced;
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

    private static long median(Set<Feature2D> pixelList, boolean is1) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (Feature2D px : pixelList) {
            if (is1) {
                stats.addValue(px.getMidPt1());
            } else {
                stats.addValue(px.getMidPt2());
            }
        }
        return (long) stats.getPercentile(50);
    }

    public static long distance(long x, long y) {
        return (long) Math.sqrt(x * x + y * y);
    }
}
