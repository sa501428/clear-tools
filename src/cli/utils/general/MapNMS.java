package cli.utils.general;

import java.util.*;

public class MapNMS {

    public static final int NUM_ROWS = 5;
    public static final int COUNT_INDEX = 0;
    public static final int OE_INDEX = 1;
    public static final int PERC_INDEX = 2;
    public static final int ZSCORE_INDEX = 3;
    public static final int TOTALS_INDEX = 4;
    public static int minLocalNeeded = 5;
    public static float withinMaxThreshold = 0.95f;

    public static void populateAfterNonMaxSuppression(float[][] upMatrix, float[][] downMatrix,
                                                      Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap, int numEntries) {

        int[] upTotals = new int[numEntries];
        int[] downTotals = new int[numEntries];

        initOERows(upMatrix, downMatrix);

        Map<Integer, Map<Integer, QuadContactRecord>> cutoffMap = generateCutoffs(contactMap);
        populateAnchorScores(upMatrix, downMatrix, upTotals, downTotals, contactMap, cutoffMap);

        inPlaceDivide(upMatrix[ZSCORE_INDEX], upTotals);
        inPlaceDivide(downMatrix[ZSCORE_INDEX], downTotals);
    }

    private static void inPlaceDivide(float[] data, int[] totals) {
        for (int i = 0; i < data.length; i++) {
            if (totals[i] > 0) {
                data[i] /= totals[i];
            } else {
                data[i] = 0;
            }
        }
    }

    private static void populateAnchorScores(float[][] upMatrix, float[][] downMatrix, int[] upTotals, int[] downTotals,
                                             Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap,
                                             Map<Integer, Map<Integer, QuadContactRecord>> cutoffMap) {
        for (int i : cutoffMap.keySet()) {
            for (int j : cutoffMap.get(i).keySet()) {
                QuadContactRecord cutoff = cutoffMap.get(i).get(j);
                for (QuadContactRecord qcr : contactMap.get(i).get(j)) {
                    if (qcr.getCounts() >= cutoff.getCounts() || qcr.getOE() >= cutoff.getOE() ||
                            qcr.getPerc() >= cutoff.getPerc() || qcr.getZscore() >= cutoff.getZscore()) {

                        fillInScoreRows(upMatrix, downMatrix,
                                qcr.getBinX(), qcr.getBinY(),
                                qcr.getCounts(), qcr.getOE(), qcr.getPerc(), qcr.getZscore());
                    }
                }
            }
        }
    }


    public static Map<Integer, Map<Integer, QuadContactRecord>> generateCutoffs(Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap) {

        Map<Integer, Map<Integer, QuadContactRecord>> cutoffMap = new HashMap<>();

        for (Integer i : contactMap.keySet()) {
            for (Integer j : contactMap.get(i).keySet()) {
                List<QuadContactRecord> list = addNeighborsToListIfPossible(i, j, contactMap);
                if (list.size() > minLocalNeeded) {
                    float[] maxValues = getMaxThresholds(list);

                    if (!cutoffMap.containsKey(i)) {
                        cutoffMap.put(i, new HashMap<>());
                    }
                    cutoffMap.get(i).put(j, new QuadContactRecord(i, j,
                            withinMaxThreshold * maxValues[COUNT_INDEX],
                            withinMaxThreshold * maxValues[OE_INDEX],
                            withinMaxThreshold * maxValues[PERC_INDEX],
                            withinMaxThreshold * maxValues[ZSCORE_INDEX]));
                }

            }
        }
        return cutoffMap;
    }

    private static float[] getMaxThresholds(List<QuadContactRecord> list) {
        float[] maxValues = new float[4];
        maxValues[COUNT_INDEX] = list.get(0).getCounts();
        maxValues[OE_INDEX] = list.get(0).getOE();
        maxValues[PERC_INDEX] = list.get(0).getPerc();
        maxValues[ZSCORE_INDEX] = list.get(0).getZscore();

        for (QuadContactRecord qcr : list) {
            if (qcr.getCounts() > maxValues[COUNT_INDEX]) {
                maxValues[COUNT_INDEX] = qcr.getCounts();
            }
            if (qcr.getOE() > maxValues[OE_INDEX]) {
                maxValues[OE_INDEX] = qcr.getOE();
            }
            if (qcr.getPerc() > maxValues[PERC_INDEX]) {
                maxValues[PERC_INDEX] = qcr.getPerc();
            }
            if (qcr.getZscore() > maxValues[ZSCORE_INDEX]) {
                maxValues[ZSCORE_INDEX] = qcr.getZscore();
            }
        }
        return maxValues;
    }

    private static List<QuadContactRecord> addNeighborsToListIfPossible(Integer i0,
                                                                        Integer j0,
                                                                        Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap) {
        List<QuadContactRecord> list = new ArrayList<>(20);
        for (int i = i0 - 1; i <= i0 + 1; i++) {
            for (int j = j0 - 1; j <= j0 + 1; j++) {
                if (contactMap.containsKey(i) && contactMap.get(i).containsKey(j)) {
                    list.addAll(contactMap.get(i).get(j));
                }
            }
        }
        return list;
    }

    public static void initOERows(float[][] up, float[][] down) {
        Arrays.fill(up[OE_INDEX], 1);
        Arrays.fill(down[OE_INDEX], 1);
    }

    public static void fillInScoreRows(float[][] upMatrix, float[][] downMatrix,
                                       int binX, int binY, float counts, float oe, float perc, float zscore) {
        fillInScoreRows2(upMatrix, binX, counts, oe, perc, zscore);
        fillInScoreRows2(downMatrix, binY, counts, oe, perc, zscore);
    }

    private static void fillInScoreRows2(float[][] matrix, int bin, float counts, float oe, float perc, float zscore) {
        matrix[COUNT_INDEX][bin] += counts;
        matrix[OE_INDEX][bin] *= oe;
        matrix[PERC_INDEX][bin] += perc;
        matrix[ZSCORE_INDEX][bin] += zscore;
        matrix[TOTALS_INDEX][bin]++;
    }
}
