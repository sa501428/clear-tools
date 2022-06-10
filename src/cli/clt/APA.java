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
import cli.utils.apa.APARegionStatistics;
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
 * <p/>
 * Implemented in Juicer by mshamim
 * <p/>
 * ---
 * APA
 * ---
 * The "apa" command takes three required arguments and a number of optional
 * arguments.
 * <p/>
 * apa [-n minval] [-x maxval] [-w window]  [-r resolution(s)] [-c chromosome(s)]
 * [-k NONE/VC/VC_SQRT/KR] <HiC file(s)> <PeaksFile> <SaveFolder> [SavePrefix]
 * <p/>
 * The required arguments are:
 * <p/>
 * <hic file(s)>: Address of HiC File(s) which should end with ".hic". This is the file you will
 * load into Juicebox. URLs or local addresses may be used. To sum multiple HiC Files together,
 * use the '+' symbol between the addresses (no whitespace between addresses)
 * <PeaksFile>: List of peaks in standard 2D feature format (chr1 x1 x2 chr2 y1 y2 color ...)
 * <SaveFolder>: Working directory where outputs will be saved
 * <p/>
 * The optional arguments are:
 * -n <int> minimum distance away from the diagonal. Used to filter peaks too close to the diagonal.
 * Units are in terms of the provided resolution. (e.g. -n 30 @ resolution 5kB will filter loops
 * within 30*(5000/sqrt(2)) units of the diagonal)
 * -x <int> maximum distance away from the diagonal. Used to filter peaks too far from the diagonal.
 * Units are in terms of the provided resolution. (e.g. -n 30 @ resolution 5kB will filter loops
 * further than 30*(5000/sqrt(2)) units of the diagonal)
 * -r <int(s)> resolution for APA; multiple resolutions can be specified using commas (e.g. 5000,10000)
 * -c <String(s)> Chromosome(s) on which APA will be run. The number/letter for the chromosome can be
 * used with or without appending the "chr" string. Multiple chromosomes can be specified using
 * commas (e.g. 1,chr2,X,chrY)
 * -k <NONE/VC/VC_SQRT/KR> Normalizations (case sensitive) that can be selected. Generally,
 * KR (Knight-Ruiz) balancing should be used when available.
 * <p/>
 * Default settings of optional arguments:
 * -n 30
 * -x (infinity)
 * -r 25000,10000
 * -c (all_chromosomes)
 * -k KR
 * <p/>
 * ------------
 * APA Examples
 * ------------
 * <p/>
 * apa HIC006.hic all_loops.txt results1
 * > This command will run APA on HIC006 using loops from the all_loops files
 * > and save them under the results1 folder.
 * <p/>
 * apa https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic
 * all_loops.txt results1
 * > This command will run APA on the GM12878 mega map using loops from the all_loops
 * > files and save them under the results1 folder.
 * <p/>
 * apa -r 10000,5000 -c 17,18 HIC006.hic+HIC007.hic all_loops.txt results
 * > This command will run APA at 50 kB resolution on chromosomes 17 and 18 for the
 * > summed HiC maps (HIC006 and HIC007) using loops from the all_loops files
 * > and save them under the results folder
 */
public class APA {
    private final Object key = new Object();
    private boolean dontIncludePlots = false;
    private String loopListPath;
    private File outputDirectory;
    private Dataset ds;
    private NormalizationType norm;
    //defaults
    // TODO right now these units are based on n*res/sqrt(2)
    // TODO the sqrt(2) scaling should be removed (i.e. handle scaling internally)
    private int minPeakDist = 30; // distance between two bins, can be changed in opts
    private int maxPeakDist = Integer.MAX_VALUE;
    private int window = 10;
    private int numCPUThreads = 4;
    private int resolution = 5000;
    private int regionWidth = 10;
    private boolean includeInterChr = false;
    private boolean aggregateNormalization = false;

    public static String getBasicUsage() {
        return "apa <hicFile(s)> <PeaksFile> <SaveFolder>";
    }

    protected void readArguments(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        loopListPath = args[2];
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, true);

        norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{parser.getNormalizationStringOption(), "SCALE", "KR", "NONE"});

        int potentialWindow = parser.getWindowSizeOption(10);
        if (potentialWindow > 0)
            window = potentialWindow;

        int potentialMinPeakDist = parser.getMinDistVal(2 * window);
        if (potentialMinPeakDist > -1)
            minPeakDist = potentialMinPeakDist;

        int potentialMaxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);
        if (potentialMaxPeakDist > 0)
            maxPeakDist = potentialMaxPeakDist;

        includeInterChr = parser.getIncludeInterChromosomal();

        int possibleRegionWidth = parser.getCornerRegionDimensionOption(window / 2);
        if (possibleRegionWidth > 0)
            regionWidth = possibleRegionWidth;

        resolution = parser.getResolutionOption(5000);

        numCPUThreads = parser.getNumThreads(4);

        aggregateNormalization = parser.getAggregateNormalization();
    }

    private void printUsageAndExit() {
        System.out.println("apa [-n minval] [-x maxval] [-w window] [-r resolution(s)] [-c chromosomes]" +
                " [-k NONE/VC/VC_SQRT/KR] [-q corner_width] [-e include_inter_chr] [-u save_all_data] [--aggnorm]" +
                " <hicFile(s)> <PeaksFile> <SaveFolder>");
        System.exit(19);
    }

    public void run() {
        APARegionStatistics result = null;
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
                                        System.err.println(e);
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
            //save data as int array
            result = APADataStack.retrieveDataStatistics(regionWidth); //should retrieve data
            Integer[] gwPeakNumbersArray = {gwPeakNumbers[0].get(), gwPeakNumbers[1].get(), gwPeakNumbers[2].get()};
            APADataStack.exportGenomeWideData(gwPeakNumbersArray, regionWidth, dontIncludePlots);
            APADataStack.clearAllData();
        } else {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }


        System.out.println("APA complete");
    }
}