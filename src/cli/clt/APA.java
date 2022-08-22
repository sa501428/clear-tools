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

package cli.clt;

import cli.utils.apa.APADataExporter;
import cli.utils.apa.APAUtils;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import cli.utils.general.HiCUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Aggregate Peak Analysis developed by mhuntley
 * Implemented in Juicer by mshamim
 * Various updates by mshamim and suhas-rao
 */
public class APA {
    private final String loopListPath;
    private final File outputDirectory;
    private final Dataset ds;
    private NormalizationType norm;
    //defaults
    // TODO right now these units are based on n*res/sqrt(2)
    // TODO the sqrt(2) scaling should be removed (i.e. handle scaling internally)
    private final int minPeakDist; // distance between two bins, can be changed in opts
    private final int maxPeakDist;
    private final int window;
    private final int matrixWidthL;
    private final int numCPUThreads;
    private final int resolution;
    private final boolean includeInterChr;
    private final boolean useAgNorm;

    private final float[][] globalAPAMatrix;
    private final double[] globalRowSum;
    private final double[] globalColSum;

    public APA(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, true);
        loopListPath = args[2];
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        useAgNorm = parser.getAggregateNormalization();
        if (useAgNorm) {
            norm = NormalizationHandler.NONE;
        } else {
            String possibleNorm = parser.getNormalizationStringOption();
            try {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } catch (Exception e) {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
            }
        }
        System.out.println("Using normalization: " + norm.getLabel());
        window = parser.getWindowSizeOption(10);
        minPeakDist = parser.getMinDistVal(2 * window);
        maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);
        includeInterChr = parser.getIncludeInterChromosomal();
        resolution = parser.getResolutionOption(5000);
        numCPUThreads = parser.getNumThreads(4);

        matrixWidthL = 2 * window + 1;
        globalAPAMatrix = new float[matrixWidthL][matrixWidthL];
        if (useAgNorm) {
            globalRowSum = new double[matrixWidthL];
            globalColSum = new double[matrixWidthL];
        } else {
            globalRowSum = null;
            globalColSum = null;
        }
    }

    private void printUsageAndExit() {
        System.out.println("apa [--min-dist minval] [--max-dist max_val] [--window window] [-r resolution]" +
                " [-k NONE/VC/VC_SQRT/KR] [--corner-width corner_width] [--include-inter include_inter_chr] [--ag-norm]" +
                " <input.hic> <loops.bedpe> <outfolder>");
        System.exit(19);
    }

    public void run() {
        System.out.println("Processing APA for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(resolution);

        ChromosomeHandler handler = ds.getChromosomeHandler();
        final AtomicInteger[] gwPeakNumbers = {new AtomicInteger(0), new AtomicInteger(0), new AtomicInteger(0)};

        Feature2DList loopList = loadLoopsAPAStyle(gwPeakNumbers, handler);
        if (loopList.getNumTotalFeatures() < 1) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        int pairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), includeInterChr);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger currNumLoops = new AtomicInteger(0);

        final NormalizationType vcNorm = ds.getNormalizationHandler().getNormTypeFromString("VC");

        final Object key = new Object();

        ParallelizationTools.launchParallelizedCode(() -> {

            float[][] output = new float[matrixWidthL][matrixWidthL];
            double[] rowSum = null;
            double[] colSum = null;

            if (useAgNorm) {
                rowSum = new double[matrixWidthL];
                colSum = new double[matrixWidthL];
            }

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < pairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                Matrix matrix = ds.getMatrix(chr1, chr2);
                if (matrix != null) {


                    List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                    if (loops != null && loops.size() > 0) {

                        double[] vector1 = null;
                        double[] vector2 = null;
                        if (useAgNorm) {
                            vector1 = ds.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);

                            if (chr1.getIndex() == chr2.getIndex()) {
                                vector2 = vector1;
                            } else {
                                vector2 = ds.getNormalizationVector(chr2.getIndex(), zoom, vcNorm).getData().getValues().get(0);
                            }
                        }

                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {
                                for (Feature2D loop : loops) {

                                    Utils.addLocalizedData(output, zd, loop, matrixWidthL, resolution, window, norm, key);
                                    if (useAgNorm) {
                                        int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
                                        int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
                                        APAUtils.addLocalRowSums(rowSum, vector1, binXStart);
                                        APAUtils.addLocalRowSums(colSum, vector2, binYStart);
                                    }

                                    if (currNumLoops.incrementAndGet() % 100 == 0) {
                                        System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                    }
                                }
                                System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                            } catch (Exception e) {
                                System.err.println(e.getMessage());
                            }
                        }
                        vector1 = null;
                        vector2 = null;
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }

            synchronized (key) {
                APAUtils.inPlaceSumMatrices(globalAPAMatrix, output);
                if (useAgNorm) {
                    APAUtils.inPlaceSumVectors(globalRowSum, rowSum);
                    APAUtils.inPlaceSumVectors(globalColSum, colSum);
                }
            }
        });

        System.out.println("Exporting APA results...");
        APADataExporter.exportGenomeWideData(gwPeakNumbers, outputDirectory, useAgNorm, globalAPAMatrix,
                globalRowSum, globalColSum);
        System.out.println("APA complete");
    }

    private Feature2DList loadLoopsAPAStyle(AtomicInteger[] gwPeakNumbers, ChromosomeHandler handler) {
        return Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> {
                    List<Feature2D> uniqueFeatures = new ArrayList<>(new HashSet<>(features));
                    List<Feature2D> filteredUniqueFeatures = APAUtils.filterFeaturesBySize(uniqueFeatures,
                            minPeakDist, maxPeakDist, resolution);
                    gwPeakNumbers[0].addAndGet(filteredUniqueFeatures.size());
                    gwPeakNumbers[1].addAndGet(uniqueFeatures.size());
                    gwPeakNumbers[2].addAndGet(features.size());
                    return filteredUniqueFeatures;
                }, false);
    }
}