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

package cli.clt.apa;

import cli.clt.CommandLineParser;
import cli.utils.apa.APADataExporter;
import cli.utils.apa.APAUtils;
import cli.utils.data.Bounds2DInfo;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class APA3 {
    public static String usage = "apa3 [--ag-norm] [-k NORM] [--window val]" +
            " [--min-dist val] [--max-dist val] [--include-inter] [-r resolution]" +
            " <input.hic> <loops1.bedpe> <output1.npy> [<loops2.bedpe> <output2.npy> ...]";
    protected final NormalizationType norm;
    protected final int window;
    protected final int matrixWidthL;
    protected final int resolution;
    private final Dataset ds;
    //defaults
    // TODO right now these units are based on n*res/sqrt(2)
    // TODO the sqrt(2) scaling should be removed (i.e. handle scaling internally)
    private final int minPeakDist; // distance between two bins, can be changed in opts
    private final int maxPeakDist;
    private final boolean includeInterChr;
    private final boolean useAgNorm;
    private final int numOutputs;

    private final String[] outputFilePaths;

    private final List<float[][]> globalAPAMatrices;
    private final List<double[]> globalRowSums;
    private final List<double[]> globalColSums;

    private final AtomicInteger currChromPair = new AtomicInteger(0);
    private final NormalizationType vcNorm = NormalizationHandler.VC;
    private final Object key = new Object();
    private final HiCZoom zoom;
    private final ChromosomeHandler handler;
    private final Feature2DList[] allLoopLists;
    private int numThreads = Runtime.getRuntime().availableProcessors();

    public APA3(String[] args, CommandLineParser parser) {
        if (args.length < 4) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(5000);
        boolean useBI = resolution >= 50;

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, useBI);
        numOutputs = (args.length - 2) / 2;
        String[] loopListPaths = new String[numOutputs];
        outputFilePaths = new String[numOutputs];
        int counter = 0;
        for (int i = 2; i < args.length - 1; i += 2) {
            loopListPaths[counter] = args[i];
            outputFilePaths[counter] = args[i + 1];
            counter++;
            // HiCFileTools.createValidDirectory(args[3]);
        }

        String possibleNorm = parser.getNormalizationStringOption();
        useAgNorm = parser.getAggregateNormalization() || isAgNorm(possibleNorm);

        NormalizationType tempNorm = NormalizationHandler.NONE;
        if (!useAgNorm) {
            try {
                tempNorm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } catch (Exception e) {
                tempNorm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
            }
        }

        numThreads = parser.getNumThreads(numThreads);

        norm = tempNorm;
        System.out.println("Using normalization: " + norm.getLabel());
        if (useAgNorm) {
            System.out.println("Will apply aggregate normalization.");
        }
        window = parser.getWindowSizeOption(10);
        minPeakDist = parser.getMinDistVal(2 * window);
        maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);
        includeInterChr = parser.getIncludeInterChromosomal();

        matrixWidthL = 2 * window + 1;
        globalAPAMatrices = initializeMatrices(numOutputs, matrixWidthL);
        globalRowSums = initializeArrays(numOutputs, matrixWidthL);
        globalColSums = initializeArrays(numOutputs, matrixWidthL);

        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();

        allLoopLists = loadLoopsAPAStyle(loopListPaths, handler);
    }

    private List<double[]> initializeArrays(int numOutputs, int matrixWidthL) {
        List<double[]> arrays = new ArrayList<>(numOutputs);
        for (int i = 0; i < numOutputs; i++) {
            arrays.add(new double[matrixWidthL]);
        }
        return arrays;
    }

    private List<float[][]> initializeMatrices(int numOutputs, int matrixWidthL) {
        List<float[][]> matrices = new ArrayList<>(numOutputs);
        for (int i = 0; i < numOutputs; i++) {
            matrices.add(new float[matrixWidthL][matrixWidthL]);
        }
        return matrices;
    }

    private boolean isAgNorm(String norm) {
        String normLower = norm.toLowerCase();
        return normLower.contains("ag") && normLower.contains("norm");
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }

    public void run() {
        System.out.println("Processing APA for resolution " + resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new HashMap<>();
        int pairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), includeInterChr);

        ParallelizationTools.launchParallelizedCode(numThreads, () -> {

            List<float[][]> outputs = initializeMatrices(numOutputs, matrixWidthL);
            List<double[]> rowSums = initializeArrays(numOutputs, matrixWidthL);
            List<double[]> colSums = initializeArrays(numOutputs, matrixWidthL);

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < pairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                Matrix matrix = ds.getMatrix(chr1, chr2);
                if (matrix != null) {
                    List<List<Feature2D>> loops = getRelevantLoopLists(allLoopLists, chr1.getIndex(), chr2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {
                                processLoopsForRegion(zd, loops, outputs);
                                if (useAgNorm) {
                                    doAggregateNormalization(chr1, chr2, zoom, vcNorm, loops, rowSums, colSums);
                                }
                                System.out.println("Completed " + chr1.getName() + "_" + chr2.getName());
                            } catch (Exception e) {
                                System.err.println(e.getMessage());
                            }
                        }
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }

            synchronized (key) {
                for (int i = 0; i < numOutputs; i++) {
                    APAUtils.inPlaceSumMatrices(globalAPAMatrices.get(i), outputs.get(i));
                    if (useAgNorm) {
                        APAUtils.inPlaceSumVectors(globalRowSums.get(i), rowSums.get(i));
                        APAUtils.inPlaceSumVectors(globalColSums.get(i), colSums.get(i));
                    }
                }
            }
        });

        System.out.println("Exporting APA results...");
        for (int i = 0; i < numOutputs; i++) {
            APADataExporter.simpleExportGenomeWideData(outputFilePaths[i], useAgNorm, globalAPAMatrices.get(i),
                    globalRowSums.get(i), globalColSums.get(i));
        }
        System.out.println("APA complete");
    }

    private List<List<Feature2D>> getRelevantLoopLists(Feature2DList[] allLoopLists, int index1, int index2) {
        List<List<Feature2D>> relevantLoops = new ArrayList<>(allLoopLists.length);
        for (Feature2DList loopList : allLoopLists) {
            relevantLoops.add(loopList.get(index1, index2));
        }
        return relevantLoops;
    }

    private void doAggregateNormalization(Chromosome chr1, Chromosome chr2, HiCZoom zoom,
                                          NormalizationType vcNorm,
                                          List<List<Feature2D>> cloops, List<double[]> rowSums, List<double[]> colSums) {
        double[] vector1 = ds.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);
        double[] vector2 = vector1;
        if (chr1.getIndex() != chr2.getIndex()) {
            vector2 = ds.getNormalizationVector(chr2.getIndex(), zoom, vcNorm).getData().getValues().get(0);
        }

        for (int i = 0; i < cloops.size(); i++) {
            List<Feature2D> loops = cloops.get(i);
            double[] rowSum = rowSums.get(i);
            double[] colSum = colSums.get(i);

            for (Feature2D loop : loops) {
                int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
                int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
                APAUtils.addLocalRowSums(rowSum, vector1, binXStart);
                APAUtils.addLocalRowSums(colSum, vector2, binYStart);
            }
        }
    }

    protected void processLoopsForRegion(MatrixZoomData zd, List<List<Feature2D>> allLoops,
                                         List<float[][]> outputs) {
        RegionsOfInterest roi = new RegionsOfInterest(resolution, window, matrixWidthL, allLoops);
        int counter = 0;

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (roi.probablyContainsRecord(cr)) {
                    for (int i = 0; i < allLoops.size(); i++) {
                        if (roi.containsRecord(cr, i)) {
                            populateMatrixIfApplicable(outputs.get(i), cr, roi.getBoundsInfo(cr, i));
                        }
                    }
                    if (counter++ % 10000 == 0) {
                        System.out.print(".");
                    }
                }
            }
        }

        roi.clear();
    }

    private void populateMatrixIfApplicable(float[][] matrix, ContactRecord cr,
                                            List<Bounds2DInfo> allBounds) {
        for (Bounds2DInfo bound : allBounds) {
            if (bound.contains(cr.getBinX(), cr.getBinY())) {
                int relativeX = cr.getBinX() - bound.getBinXStart();
                int relativeY = cr.getBinY() - bound.getBinYStart();
                matrix[relativeX][relativeY] += cr.getCounts();
            }
        }
    }


    private Feature2DList[] loadLoopsAPAStyle(String[] loopListPath, ChromosomeHandler handler) {
        Feature2DList[] loopLists = new Feature2DList[loopListPath.length];
        for (int i = 0; i < loopListPath.length; i++) {
            loopLists[i] = loadLoopsAPAStyle(loopListPath[i], handler);
            if (loopLists[i] == null || loopLists[i].getNumTotalFeatures() < 1) {
                System.err.println("Loop list: " + loopListPath[i] + " is empty or incorrect path provided.");
                System.exit(3);
            }
        }
        return loopLists;
    }

    private Feature2DList loadLoopsAPAStyle(String loopListPath, ChromosomeHandler handler) {
        return Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> APAUtils.filterFeaturesBySize(new ArrayList<>(new HashSet<>(features)),
                        minPeakDist, maxPeakDist, resolution), false);
    }
}