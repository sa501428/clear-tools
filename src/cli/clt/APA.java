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

import cli.Main;
import cli.utils.apa.APADataStack;
import cli.utils.apa.APAUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
    private final int numCPUThreads;
    private final int resolution;
    private final boolean includeInterChr;
    private final boolean aggregateNormalization;

    public APA(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, true);
        loopListPath = args[2];
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
        }
        System.out.println("Using normalization: " + norm.getLabel());
        window = parser.getWindowSizeOption(10);
        minPeakDist = parser.getMinDistVal(2 * window);
        maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);
        includeInterChr = parser.getIncludeInterChromosomal();
        resolution = parser.getResolutionOption(5000);
        numCPUThreads = parser.getNumThreads(4);
        aggregateNormalization = parser.getAggregateNormalization();
    }

    private void printUsageAndExit() {
        System.out.println("apa [--min-dist minval] [--max-dist max_val] [--window window] [-r resolution]" +
                " [-k NONE/VC/VC_SQRT/KR] [--corner-width corner_width] [--include-inter include_inter_chr] [--ag-norm]" +
                " <input.hic> <loops.bedpe> <outfolder>");
        System.exit(19);
    }

    public void run() {
        int L = 2 * window + 1;

        AtomicInteger[] gwPeakNumbers = {new AtomicInteger(0), new AtomicInteger(0), new AtomicInteger(0)};

        System.out.println("Processing APA for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(resolution);

        ChromosomeHandler handler = ds.getChromosomeHandler();

        // Metrics resulting from apa filtering
        final Map<String, Integer[]> filterMetrics = new HashMap<>();
        //looplist is empty here why??
        // Remove duplicates and filters by size
        // also save internal metrics for these measures
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> {
                    List<Feature2D> uniqueFeatures = new ArrayList<>(new HashSet<>(features));
                    List<Feature2D> filteredUniqueFeatures = APAUtils.filterFeaturesBySize(uniqueFeatures,
                            minPeakDist, maxPeakDist, resolution);

                    filterMetrics.put(chr,
                            new Integer[]{filteredUniqueFeatures.size(), uniqueFeatures.size(), features.size()});

                    return filteredUniqueFeatures;
                }, false);

        if (loopList.getNumTotalFeatures() > 0) {

            double maxProgressStatus = handler.size();
            final AtomicInteger currentProgressStatus = new AtomicInteger(0);
            List<Chromosome[]> chromosomePairs = new ArrayList<>();
            for (Chromosome chr1 : handler.getChromosomeArrayWithoutAllByAll()) {
                for (Chromosome chr2 : handler.getChromosomeArrayWithoutAllByAll()) {
                    Chromosome[] chromosomePair = {chr1, chr2};
                    chromosomePairs.add(chromosomePair);
                }
            }

            APADataStack.initializeDataSaveFolder(outputDirectory, "" + resolution);

            for (int chrPair = 0; chrPair < chromosomePairs.size(); chrPair++) {
                Chromosome[] pair = chromosomePairs.get(chrPair);
                Chromosome chr1 = pair[0];
                Chromosome chr2 = pair[1];
                if ((chr1.getIndex() < chr2.getIndex() && includeInterChr) || (chr2.getIndex() == chr1.getIndex())) {
                    MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, zoom);
                    if (zd == null) continue;
                    List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                    if (loops == null || loops.size() == 0) {
                        if (Main.printVerboseComments) {
                            System.out.println("CHR " + chr1.getName() + " CHR " + chr2.getName() + " - no loops, check loop filtering constraints");
                        }
                        continue;
                    } else if (Main.printVerboseComments) {
                        System.out.println("CHR " + chr1.getName() + " " + chr1.getIndex() + " CHR " + chr2.getName() + " " + chr2.getIndex());
                    }

                    int numOfLoopChunks = (loops.size() / 2) + 1;
                    int numOfLoops = loops.size();
                    final AtomicInteger loopChunk = new AtomicInteger(0);
                    Integer[] peakNumbers = filterMetrics.get(Feature2DList.getKey(chr1, chr2));

                    if (loops.size() != peakNumbers[0]) {
                        System.err.println("Error reading statistics from " + chr1 + chr2);
                    }

                    for (int i = 0; i < peakNumbers.length; i++) {
                        gwPeakNumbers[i].addAndGet(peakNumbers[i]);
                    }

                    ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
                    for (int l = 0; l < numCPUThreads; l++) {
                        final int threadChrPair = chrPair;
                        Runnable worker = () -> {
                            APADataStack apaDataStack = new APADataStack(L, chromosomePairs.size(), aggregateNormalization);
                            int threadChunk = loopChunk.getAndIncrement();
                            while (threadChunk < numOfLoopChunks) {
                                for (int loopIndex = threadChunk * 2; loopIndex < Math.min(numOfLoops, (threadChunk + 1) * 2); loopIndex++) {
                                    Feature2D loop = loops.get(loopIndex);
                                    try {
                                        RealMatrix newData = APAUtils.extractLocalizedData(zd, loop, L, resolution, window, norm);
                                        apaDataStack.addData(newData);
                                        if (aggregateNormalization) {
                                            List<RealMatrix> newVectors = APAUtils.extractLocalizedRowSums(zd, loop, L, resolution, window, norm);
                                            apaDataStack.addRowSums(newVectors);
                                        }
                                    } catch (Exception e) {
                                        System.err.println(e.getLocalizedMessage());
                                        System.err.println("Unable to find data for loop: " + loop);
                                    }
                                }
                                threadChunk = loopChunk.getAndIncrement();
                            }
                            synchronized (ds) {
                                apaDataStack.updateGenomeWideData();
                            }

                            apaDataStack.updateChromosomeWideData(threadChrPair);
                        };
                        executor.execute(worker);
                    }
                    executor.shutdown();

                    // Wait until all threads finish
                    while (!executor.isTerminated()) {
                    }
                }
                if (chr2.getIndex() == chr1.getIndex()) {
                    System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus)) + "% ");
                }
            }

            System.out.println("Exporting APA results...");
            Integer[] gwPeakNumbersArray = {gwPeakNumbers[0].get(), gwPeakNumbers[1].get(), gwPeakNumbers[2].get()};
            APADataStack.exportGenomeWideData(gwPeakNumbersArray);
            APADataStack.clearAllData();
        } else {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }


        System.out.println("APA complete");
    }
}