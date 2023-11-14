package cli.utils.anchors;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class LoopToAnchorIntersecter {
    public static Feature2DList intersect(List<Feature2DList> loopLists, ChromosomeHandler handler,
                                          Map<String, float[]> chromToHiResAnchorCounts, int highResolution) {

        Map<String, Map<Integer, Set<Integer>>> globalAnswer = new HashMap<>();
        Map<String, String[]> keyToChroms = new HashMap<>();

        int hiResWindow = 5;

        AtomicInteger globalIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {

            Map<String, Map<Integer, Set<Integer>>> localAnswer = new HashMap<>();

            int i = globalIndex.getAndIncrement();
            while (i < loopLists.size()) {
                Feature2DList loopList = loopLists.get(i);

                loopList.processLists((sKey, list) -> {

                    synchronized (keyToChroms) {
                        if (list.size() > 0) {
                            keyToChroms.put(sKey, new String[]{list.get(0).getChr1(), list.get(0).getChr2()});
                        }
                    }

                    for (Feature2D loop : list) {

                        int x = getIdealPosition(loop.getAttribute("localX"), highResolution,
                                chromToHiResAnchorCounts.get(loop.getChr1()), hiResWindow,
                                loop.getStart1(), loop.getEnd1());
                        int y = getIdealPosition(loop.getAttribute("localY"), highResolution,
                                chromToHiResAnchorCounts.get(loop.getChr2()), hiResWindow,
                                loop.getStart2(), loop.getEnd2());

                        if (x > 0 && y > 0) {
                            if (!localAnswer.containsKey(sKey)) {
                                localAnswer.put(sKey, new HashMap<>());
                            }
                            if (!localAnswer.get(sKey).containsKey(x)) {
                                localAnswer.get(sKey).put(x, new HashSet<>());
                            }
                            localAnswer.get(sKey).get(x).add(y);
                        }
                    }
                });

                i = globalIndex.getAndIncrement();
            }

            synchronized (globalAnswer) {
                for (String chr1 : localAnswer.keySet()) {
                    if (!globalAnswer.containsKey(chr1)) {
                        globalAnswer.put(chr1, new HashMap<>());
                    }
                    for (int x : localAnswer.get(chr1).keySet()) {
                        if (!globalAnswer.get(chr1).containsKey(x)) {
                            globalAnswer.get(chr1).put(x, new HashSet<>());
                        }
                        globalAnswer.get(chr1).get(x).addAll(localAnswer.get(chr1).get(x));
                    }
                }
            }

        });

        return globalAnswerToFeature2DList(keyToChroms, highResolution, globalAnswer);
    }

    private static Feature2DList globalAnswerToFeature2DList(Map<String, String[]> keyToChroms, int resolution,
                                                             Map<String, Map<Integer, Set<Integer>>> globalAnswer) {
        Feature2DList result = new Feature2DList();
        for (String sKey : globalAnswer.keySet()) {
            String chr1 = keyToChroms.get(sKey)[0];
            String chr2 = keyToChroms.get(sKey)[1];

            List<Feature2D> loops = new LinkedList<>();
            for (int x : globalAnswer.get(sKey).keySet()) {
                for (int y : globalAnswer.get(sKey).get(x)) {
                    long start1 = (long) x * resolution;
                    long end1 = start1 + resolution;
                    long start2 = (long) y * resolution;
                    long end2 = start2 + resolution;

                    loops.add(new Feature2D(Feature2D.FeatureType.PEAK,
                            chr1, start1, end1, chr2, start2, end2, Color.BLUE, new HashMap<>()));
                }
            }
            result.addByKey(sKey, loops);
        }
        return result;
    }

    private static int getIdealPosition(String localVal, int highResolution, float[] data, int hiResWindow,
                                        long startLowRes, long endLowRes) {

        try {
            int pos = (int) (Long.parseLong(localVal) / highResolution);
            int startCheckPos = Math.max(1, pos - hiResWindow);
            int endCheckPos = Math.min(data.length - 1, pos + hiResWindow + 1);
            pos = getBestPositionInRange(startCheckPos, endCheckPos, data, pos);

            if (isSimpleLocalMax(pos, data)) {
                return pos;
            }

        } catch (Exception ignored) {
        }

        int startCheckPos = (int) (startLowRes / highResolution);
        int endCheckPos = (int) (endLowRes / highResolution);
        startCheckPos = Math.max(1, startCheckPos - (2 * hiResWindow));
        endCheckPos = Math.min(data.length - 1, endCheckPos + (2 * hiResWindow) + 1);

        int pos = getBestPositionInRange(startCheckPos, endCheckPos, data, (startCheckPos + endCheckPos) / 2);
        if (isSimpleLocalMax(pos, data)) {
            return pos;
        }
        return -1;
    }

    private static boolean isSimpleLocalMax(int pos, float[] data) {
        return data[pos] > 1 && data[pos] > data[pos - 1] && data[pos] > data[pos + 1];
    }

    private static int getBestPositionInRange(int startCheckPos, int endCheckPos, float[] data, int initialPos) {
        int pos = initialPos;
        for (int x = startCheckPos; x < endCheckPos; x++) {
            if (data[x] > data[pos]) {
                pos = x;
            }
        }
        return pos;
    }
}
