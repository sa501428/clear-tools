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
import cli.utils.general.HiCUtils;
import cli.utils.general.Utils;
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
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class APA {
    public static String usage = "apa[2] [--ag-norm] [-k NORM] [--window val]" +
            " [--min-dist val] [--max-dist val] [--include-inter] [-r resolution]" +
            " <input.hic> <loops.bedpe> <outfolder>";
    private final String loopListPath;
    private final File outputDirectory;
    private final Dataset ds;
    protected final NormalizationType norm;
    //defaults
    // TODO right now these units are based on n*res/sqrt(2)
    // TODO the sqrt(2) scaling should be removed (i.e. handle scaling internally)
    private final int minPeakDist; // distance between two bins, can be changed in opts
    private final int maxPeakDist;
    protected final int window;
    protected final int matrixWidthL;
    protected final int resolution;
    private final boolean includeInterChr;
    private final boolean useAgNorm;

    private final float[][] globalAPAMatrix;
    private final double[] globalRowSum;
    private final double[] globalColSum;

    private final AtomicInteger currChromPair = new AtomicInteger(0);
    private final AtomicInteger currNumLoops = new AtomicInteger(0);
    private final NormalizationType vcNorm = NormalizationHandler.VC;
    private final Object key = new Object();
    private final AtomicInteger[] gwPeakNumbers = {new AtomicInteger(0),
            new AtomicInteger(0), new AtomicInteger(0)};
    private final HiCZoom zoom;
    private final ChromosomeHandler handler;
    private final Feature2DList loopList;
    private final int numTotalLoops;

    public APA(String[] args, CommandLineParser parser, boolean loadAllBlockIndices) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(5000);
        boolean useBI = !loadAllBlockIndices || resolution >= 50;

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, useBI);
        loopListPath = args[2];
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

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
        globalAPAMatrix = new float[matrixWidthL][matrixWidthL];
        globalRowSum = new double[matrixWidthL];
        globalColSum = new double[matrixWidthL];

        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();

        loopList = loadLoopsAPAStyle(gwPeakNumbers, handler);
        numTotalLoops = loopList.getNumTotalFeatures();
        if (numTotalLoops < 1) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }
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

        ParallelizationTools.launchParallelizedCode(() -> {

            float[][] output = new float[matrixWidthL][matrixWidthL];
            double[] rowSum = new double[matrixWidthL];
            double[] colSum = new double[matrixWidthL];

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < pairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                Matrix matrix = ds.getMatrix(chr1, chr2);
                if (matrix != null) {
                    List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {
                                processLoopsForRegion(zd, loops, output, currNumLoops, numTotalLoops);
                                if (useAgNorm) {
                                    doAggregateNormalization(chr1, chr2, zoom, vcNorm, loops, rowSum, colSum);
                                }
                                System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
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

    private void doAggregateNormalization(Chromosome chr1, Chromosome chr2, HiCZoom zoom, NormalizationType vcNorm, List<Feature2D> loops, double[] rowSum, double[] colSum) {
        double[] vector1 = ds.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);
        double[] vector2 = vector1;
        if (chr1.getIndex() != chr2.getIndex()) {
            vector2 = ds.getNormalizationVector(chr2.getIndex(), zoom, vcNorm).getData().getValues().get(0);
        }

        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            APAUtils.addLocalRowSums(rowSum, vector1, binXStart);
            APAUtils.addLocalRowSums(colSum, vector2, binYStart);
        }
    }

    protected void processLoopsForRegion(MatrixZoomData zd, List<Feature2D> loops, float[][] output, AtomicInteger currNumLoops, int numTotalLoops) {
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            Utils.addLocalBoundedRegion(output, zd, binXStart, binYStart, matrixWidthL, norm);
            if (currNumLoops.incrementAndGet() % 100 == 0) {
                System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
            }
        }
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