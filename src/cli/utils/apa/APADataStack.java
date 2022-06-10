/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package cli.utils.apa;


import cli.Main;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by muhammadsaadshamim on 5/1/15.
 */
public class APADataStack {

    private static final Object key = new Object();
    // genome wide variables
    private static boolean genomeWideVariablesNotSet = true;
    private static RealMatrix gwAPAMatrix;
    private static RealMatrix gwNormedAPAMatrix;
    private static RealMatrix gwCenterNormedAPAMatrix;
    private static RealMatrix gwRankAPAMatrix;
    private static List<Double> gwEnhancement;
    private static RealMatrix gwUpstreamAnchorRowSums;
    private static RealMatrix gwDownstreamAnchorRowSums;
    private static int matrixSize;
    // chromosome wide variables
    private static boolean chromosomeWideVariablesNotSet = true;
    private static Map<Integer, RealMatrix> chrAPAMatrices = new HashMap<>();
    private static Map<Integer, RealMatrix> chrNormedAPAMatrices = new HashMap<>();
    private static Map<Integer, RealMatrix> chrCenterNormedAPAMatrices = new HashMap<>();
    private static Map<Integer, RealMatrix> chrRankAPAMatrices = new HashMap<>();
    private static Map<Integer, List<Double>> chrEnhancements = new ConcurrentHashMap<>();
    private static Map<Integer, RealMatrix> chrUpstreamAnchorRowSums = new HashMap<>();
    private static Map<Integer, RealMatrix> chrDownstreamAnchorRowSums = new HashMap<>();

    // aggregate normalization variable
    private static boolean aggregateNormVariablesNotSet = true;
    private static boolean aggregateNormalization;

    // saving data variables
    private static int[] axesRange;
    private static File dataDirectory;
    private static RealMatrix upstreamAnchorRowSums;
    private static RealMatrix downstreamAnchorRowSums;
    // data stack variables
    private final List<Double> enhancement;
    private RealMatrix APAMatrix;
    private RealMatrix normedAPAMatrix;
    private RealMatrix centerNormedAPAMatrix;
    private RealMatrix rankAPAMatrix;

    /**
     * class for saving data from chromosme wide run of APA, keeps static class to store genomide data
     *
     * @param n width of matrix
     */
    public APADataStack(int n, int nChr, boolean aggNorm) {
        APAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        normedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        centerNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        rankAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        enhancement = new ArrayList<>();
        upstreamAnchorRowSums = MatrixTools.cleanArray2DMatrix(n, 1);
        downstreamAnchorRowSums = MatrixTools.cleanArray2DMatrix(n, 1);

        initializeGenomeWideVariables(n);
        initializeChromosomeWideVariables(n, nChr);
        initializeAggregateNormalizationVariables(aggNorm);
        axesRange = new int[]{-n / 2, 1, -n / 2, 1};
    }

    /**
     * Ensure that directory for saving exists
     *
     * @param outputFolderDirectory to directory
     * @param prefix                of files to be saved
     */
    public static void initializeDataSaveFolder(File outputFolderDirectory, String prefix) {
        if (prefix.length() < 1) {// no preference specified
            dataDirectory = new File(outputFolderDirectory,
                    new SimpleDateFormat("yyyy.MM.dd.HH.mm").format(new Date()));
        } else {
            dataDirectory = new File(outputFolderDirectory, prefix);
        }
        dataDirectory = HiCFileTools.createValidDirectory(dataDirectory.getAbsolutePath());
    }

    private static void initializeGenomeWideVariables(int n) {
        if (genomeWideVariablesNotSet) {
            gwAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwCenterNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwRankAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            //gwCoverage = APAUtils.cleanArray2DMatrix(n, n);
            gwEnhancement = Collections.synchronizedList(new ArrayList<>());
            gwUpstreamAnchorRowSums = MatrixTools.cleanArray2DMatrix(n, 1);
            gwDownstreamAnchorRowSums = MatrixTools.cleanArray2DMatrix(n, 1);
            matrixSize = n;
            genomeWideVariablesNotSet = false;
        }
    }

    private static void initializeChromosomeWideVariables(int n, int numOfChrPairs) {
        if (chromosomeWideVariablesNotSet) {
            for (int i = 0; i < numOfChrPairs; i++) {
                chrAPAMatrices.put(i, MatrixTools.cleanArray2DMatrix(n, n));
                chrNormedAPAMatrices.put(i, MatrixTools.cleanArray2DMatrix(n, n));
                chrCenterNormedAPAMatrices.put(i, MatrixTools.cleanArray2DMatrix(n, n));
                chrRankAPAMatrices.put(i, MatrixTools.cleanArray2DMatrix(n, n));
                chrEnhancements.put(i, Collections.synchronizedList(new ArrayList<>()));
                chrUpstreamAnchorRowSums.put(i, MatrixTools.cleanArray2DMatrix(n, 1));
                chrDownstreamAnchorRowSums.put(i, MatrixTools.cleanArray2DMatrix(n, 1));
            }
            chromosomeWideVariablesNotSet = false;
        }
    }

    private static void initializeAggregateNormalizationVariables(boolean aggNorm) {
        if (aggregateNormVariablesNotSet) {
            aggregateNormalization = aggNorm;
            aggregateNormVariablesNotSet = false;
        }
    }

    public static void exportGenomeWideData(Integer[] peakNumbers) {
        double gwNPeaksUsedInv = 1. / peakNumbers[0];
        gwNormedAPAMatrix = gwNormedAPAMatrix.scalarMultiply(gwNPeaksUsedInv);
        gwCenterNormedAPAMatrix = gwCenterNormedAPAMatrix.scalarMultiply(gwNPeaksUsedInv);
        gwRankAPAMatrix = gwRankAPAMatrix.scalarMultiply(gwNPeaksUsedInv);

        if (aggregateNormalization) {
            double gwUpstreamAnchorRowAverage = MatrixTools.getAverage(gwUpstreamAnchorRowSums);
            double gwDownstreamAnchorRowAverage = MatrixTools.getAverage(gwDownstreamAnchorRowSums);

            RealMatrix gwUpstreamAnchorRowSumsNormed = gwUpstreamAnchorRowSums.scalarMultiply(1 / gwUpstreamAnchorRowAverage);
            RealMatrix gwDownstreamAnchorRowSumsNormed = gwDownstreamAnchorRowSums.scalarMultiply(1 / gwDownstreamAnchorRowAverage);

            RealMatrix gwAggNormAPAMatrix = MatrixTools.normedCopy(gwAPAMatrix, gwUpstreamAnchorRowSumsNormed, gwDownstreamAnchorRowSumsNormed, matrixSize);
            RealMatrix[] matrices = {gwAPAMatrix, gwNormedAPAMatrix, gwCenterNormedAPAMatrix, gwRankAPAMatrix, gwAggNormAPAMatrix};
            String[] titles = {"APA", "normedAPA", "centerNormedAPA", "rankAPA", "aggNormAPA"};

            saveDataSet("gw", matrices, titles, peakNumbers);
        } else {
            RealMatrix[] matrices = {gwAPAMatrix, gwNormedAPAMatrix, gwCenterNormedAPAMatrix, gwRankAPAMatrix};
            String[] titles = {"APA", "normedAPA", "centerNormedAPA", "rankAPA"};

            saveDataSet("gw", matrices, titles, peakNumbers);
        }
    }

    public static APARegionStatistics retrieveDataStatistics(int currentRegionWidth) {
        return new APARegionStatistics(gwAPAMatrix, currentRegionWidth);
    }

    private static void saveDataSet(String prefix, RealMatrix[] apaMatrices, String[] apaDataTitles, Integer[] peakNumbers) {

        File subFolder = HiCFileTools.createValidDirectory(new File(dataDirectory, prefix).getAbsolutePath());
        if (Main.printVerboseComments) {
            System.out.println("Saving chr " + prefix + " data to " + subFolder);
        }

        String title = "N=" + peakNumbers[0] + "_(filtered)_" + peakNumbers[1] + "_(unique)_" +
                peakNumbers[2] + "_(total)";
        MatrixTools.saveMatrixTextNumpy((new File(subFolder, apaDataTitles[0] + "_" + title + ".txt")).getAbsolutePath(),
                apaMatrices[0].getData());

        if (aggregateNormalization) {
            title = "N=" + peakNumbers[0] + "_(filtered)_" + peakNumbers[1] + "_(unique)_" +
                    peakNumbers[2] + "_(total)";
            MatrixTools.saveMatrixTextNumpy((new File(subFolder, apaDataTitles[apaMatrices.length - 1] + "_" + title + ".txt")).getAbsolutePath(),
                    apaMatrices[apaMatrices.length - 1].getData());
        }
    }

    public static void clearAllData() {
        axesRange = null;
        dataDirectory = null;
        genomeWideVariablesNotSet = true;
        gwAPAMatrix = null;
        gwNormedAPAMatrix = null;
        gwCenterNormedAPAMatrix = null;
        gwRankAPAMatrix = null;
        gwEnhancement = null;
    }

    public void addData(RealMatrix newData) {
        MatrixTools.cleanUpNaNs(newData);
        APAMatrix = APAMatrix.add(newData);
        normedAPAMatrix = normedAPAMatrix.add(APAUtils.standardNormalization(newData));
        centerNormedAPAMatrix = centerNormedAPAMatrix.add(APAUtils.centerNormalization(newData));
        rankAPAMatrix = rankAPAMatrix.add(APAUtils.rankPercentile(newData));
        enhancement.add(APAUtils.peakEnhancement(newData));
    }

    public void addRowSums(List<RealMatrix> newVectors) {
        MatrixTools.cleanUpNaNs(newVectors.get(0));
        MatrixTools.cleanUpNaNs(newVectors.get(1));
        upstreamAnchorRowSums = upstreamAnchorRowSums.add(newVectors.get(0));
        downstreamAnchorRowSums = downstreamAnchorRowSums.add(newVectors.get(1));
        //System.out.println(MatrixTools.sum(upstreamAnchorRowSums.getData()) + " " + MatrixTools.sum(downstreamAnchorRowSums.getData()) + " " + MatrixTools.sum(newVectors.get(0).getData()) + " " + MatrixTools.sum(newVectors.get(1).getData()));
        //System.out.println(upstreamAnchorRowSums.transpose().getNorm() + " " + downstreamAnchorRowSums.transpose().getNorm());
    }

    public void updateGenomeWideData() {
        synchronized (key) {
            gwAPAMatrix = gwAPAMatrix.add(APAMatrix);
            gwNormedAPAMatrix = gwNormedAPAMatrix.add(normedAPAMatrix);
            gwCenterNormedAPAMatrix = gwCenterNormedAPAMatrix.add(centerNormedAPAMatrix);
            gwRankAPAMatrix = gwRankAPAMatrix.add(rankAPAMatrix);
            gwUpstreamAnchorRowSums = gwUpstreamAnchorRowSums.add(upstreamAnchorRowSums);
            gwDownstreamAnchorRowSums = gwDownstreamAnchorRowSums.add(downstreamAnchorRowSums);
        }
        synchronized (gwEnhancement) {
            gwEnhancement.addAll(enhancement);
        }
    }

    public void updateChromosomeWideData(int chrPair) {
        synchronized (chrAPAMatrices) {
            chrAPAMatrices.put(chrPair, chrAPAMatrices.get(chrPair).add(APAMatrix));
            chrNormedAPAMatrices.put(chrPair, chrNormedAPAMatrices.get(chrPair).add(normedAPAMatrix));
            chrCenterNormedAPAMatrices.put(chrPair, chrCenterNormedAPAMatrices.get(chrPair).add(centerNormedAPAMatrix));
            chrRankAPAMatrices.put(chrPair, chrRankAPAMatrices.get(chrPair).add(rankAPAMatrix));
            chrUpstreamAnchorRowSums.put(chrPair, chrUpstreamAnchorRowSums.get(chrPair).add(upstreamAnchorRowSums));
            chrDownstreamAnchorRowSums.put(chrPair, chrDownstreamAnchorRowSums.get(chrPair).add(downstreamAnchorRowSums));
        }
        synchronized (chrEnhancements) {
            chrEnhancements.get(chrPair).addAll(enhancement);
        }
    }

    public void thresholdPlots(int val) {
        MatrixTools.thresholdValues(APAMatrix, val);
    }
}
