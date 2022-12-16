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
import cli.utils.FeatureStats;
import cli.utils.apa.APAUtils;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.QuickGrouping;
import cli.utils.general.Utils;
import cli.utils.seer.SeerUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class APA1D {
    public static String usage = "apa-1d [-k NORM] [--window val]" +
            " [--min-dist val] [--max-dist val] [--include-inter] [-r resolution]" +
            " <input.hic> <loops.bedpe> <output.bedgraph>";
    private final String loopListPath;
    private final String outputPath;
    private final Dataset ds;
    //defaults
    // TODO right now these units are based on n*res/sqrt(2)
    // TODO the sqrt(2) scaling should be removed (i.e. handle scaling internally)
    private final int minPeakDist; // distance between two bins, can be changed in opts
    private final int maxPeakDist;
    private final int window;
    private final int resolution;
    private final boolean includeInterChr;
    private NormalizationType norm = NormalizationHandler.NONE;

    public APA1D(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, true);
        loopListPath = args[2];
        outputPath = args[3];

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        System.out.println("Using normalization: " + norm.getLabel());

        window = parser.getWindowSizeOption(100);
        minPeakDist = parser.getMinDistVal(3 * window);
        maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);
        includeInterChr = parser.getIncludeInterChromosomal();
        resolution = parser.getResolutionOption(10);
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }

    public void run() {
        System.out.println("Processing APA-to-1D for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(resolution);

        ChromosomeHandler handler = ds.getChromosomeHandler();

        Feature2DList loopList = loadLoopsAPAStyle(handler);
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

        // todo use double[] when not using NONE??
        final Map<Chromosome, int[]> rowSums = new HashMap<>();

        ParallelizationTools.launchParallelizedCode(() -> {

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
                                int numEntries = (int) ((chr1.getLength() / resolution) + 1);
                                int[] sums = new int[numEntries];

                                Collection<List<Feature2D>> loopGroups = QuickGrouping.groupNearbyRecords(
                                        loops, 100 * resolution).values();

                                for (List<Feature2D> group : loopGroups) {

                                    int minR = (int) ((FeatureStats.minStart1(group) / resolution) - window);
                                    int minC = (int) ((FeatureStats.minStart2(group) / resolution) - window);
                                    int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + window);
                                    int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + window);

                                    Set<ContactRecord> recordsToSum = new HashSet<>();

                                    List<ContactRecord> records = Utils.getRecords(zd, minR, minC, maxR, maxC, norm);

                                    for (Feature2D loop : group) {
                                        for (ContactRecord record : records) {
                                            if (isNearby(loop, record, resolution, window)) {
                                                recordsToSum.add(record);
                                            }
                                        }
                                    }

                                    records.clear();
                                    updateSums(sums, recordsToSum);
                                    recordsToSum.clear();

                                    if (currNumLoops.addAndGet(group.size()) % 100 == 0) {
                                        System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                    }
                                }
                                System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                                synchronized (rowSums) {
                                    rowSums.put(chr1, sums);
                                }
                            } catch (Exception e) {
                                System.err.println(e.getMessage());
                            }
                        }
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });

        System.out.println("Exporting APA results...");
        try {
            SeerUtils.exportRowSumsToBedgraph(rowSums, outputPath, resolution);
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("APA complete");
    }

    private boolean isNearby(Feature2D loop, ContactRecord record, int resolution, int window) {
        int mid1 = (int) (loop.getMidPt1() / resolution);
        int mid2 = (int) (loop.getMidPt2() / resolution);
        int dist1 = Math.abs(mid1 - record.getBinX()) + Math.abs(mid2 - record.getBinY());
        int dist2 = Math.abs(mid2 - record.getBinX()) + Math.abs(mid1 - record.getBinY());
        return dist1 < window || dist2 < window;
    }

    private void updateSums(int[] sums, Set<ContactRecord> records) {
        for (ContactRecord record : records) {
            sums[record.getBinX()] += record.getCounts();
            if (record.getBinX() != record.getBinY()) {
                sums[record.getBinY()] += record.getCounts();
            }
        }
    }

    private Feature2DList loadLoopsAPAStyle(ChromosomeHandler handler) {
        return Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> {
                    List<Feature2D> uniqueFeatures = new ArrayList<>(new HashSet<>(features));
                    return APAUtils.filterFeaturesBySize(uniqueFeatures,
                            minPeakDist, maxPeakDist, resolution);
                }, false);
    }
}