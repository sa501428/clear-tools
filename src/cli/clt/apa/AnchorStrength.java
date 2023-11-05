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
import cli.utils.general.ArrayTools;
import cli.utils.general.VectorCleaner;
import cli.utils.seer.SeerUtils;
import javastraw.expected.ExpectedUtils;
import javastraw.expected.LogExpectedZscoreSpline;
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

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorStrength {
    public static String usage = "anchor-strength [-k NORM] [-c chrom]" +
            "[--min-dist val] [--max-dist val] [-r resolution] <input.hic> <output.stem>\n" +
            "calculate localized row sums near loops";
    private final String outputPath;
    private final Dataset ds;
    private final int minPeakDist, maxPeakDist; // distance between two bins, can be changed in opts
    private final int resolution;
    private NormalizationType norm = NormalizationHandler.VC;
    private String chrom = null;

    public AnchorStrength(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(2000);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 50);
        outputPath = args[2];
        chrom = parser.getChromosomeOption();

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        System.out.println("Using normalization: " + norm.getLabel());
        minPeakDist = parser.getMinDistVal(Math.max(5, 1000 / resolution));
        maxPeakDist = parser.getMaxDistVal(10000000 / resolution);
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }

    public void run() {
        System.out.println("Processing 1D anchor strength for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(resolution);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        Chromosome[] chromosomes = getChromosomes(handler);
        final AtomicInteger currChromPair = new AtomicInteger(0);

        final Map<Chromosome, float[]> allUpStreamRowSums = new HashMap<>();
        final Map<Chromosome, float[]> allDownStreamRowSums = new HashMap<>();
        final Map<Chromosome, float[]> allUpStreamZscores = new HashMap<>();
        final Map<Chromosome, float[]> allDownStreamZscores = new HashMap<>();
        final Map<Chromosome, int[]> allUpStreamDists = new HashMap<>();
        final Map<Chromosome, int[]> allDownStreamDists = new HashMap<>();
        final Map<Chromosome, float[]> allUpStreamVals = new HashMap<>();
        final Map<Chromosome, float[]> allDownStreamVals = new HashMap<>();
        final Map<Chromosome, int[]> allUpStreamCounts = new HashMap<>();
        final Map<Chromosome, int[]> allDownStreamCounts = new HashMap<>();


        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomes.length) {
                Chromosome chrom = chromosomes[threadPair];
                Matrix matrix = ds.getMatrix(chrom, chrom);
                if (matrix != null) {

                    MatrixZoomData zd = matrix.getZoomData(zoom);
                    if (zd != null) {
                        try {
                            int numEntries = (int) ((chrom.getLength() / resolution) + 1);
                            float[] upStreamSums = new float[numEntries];
                            float[] downStreamSums = new float[numEntries];
                            int[] upStreamDists = new int[numEntries];
                            int[] downStreamDists = new int[numEntries];
                            float[] upStreamZscores = new float[numEntries];
                            float[] downStreamZscores = new float[numEntries];
                            float[] upStreamVals = new float[numEntries];
                            float[] downStreamVals = new float[numEntries];
                            int[] upStreamCounts = new int[numEntries];
                            int[] downStreamCounts = new int[numEntries];

                            LogExpectedZscoreSpline poly = new LogExpectedZscoreSpline(zd, norm, chrom, resolution);

                            double[] vector = ArrayTools.copy(ds.getNormalizationVector(chrom.getIndex(), zoom,
                                    NormalizationHandler.VC).getData().getValues().get(0));
                            VectorCleaner.inPlaceZscore(vector);

                            Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
                            while (it.hasNext()) {
                                ContactRecord cr = it.next();
                                if (cr.getCounts() > 0) {
                                    if (vector[cr.getBinX()] > -1 && vector[cr.getBinY()] > -1) {
                                        int dist = ExpectedUtils.getDist(cr);
                                        if (dist > minPeakDist && dist < maxPeakDist) {
                                            float oe = (float) ((cr.getCounts() + 1) / (poly.getExpectedFromUncompressedBin(dist) + 1));
                                            float zscore = (float) poly.getZscoreForObservedUncompressedBin(dist, cr.getCounts());
                                            if (zscore > 1 && oe > 2) {
                                                /*
                                                upStreamSums[cr.getBinX()] += (dist * zscore);
                                                upStreamDists[cr.getBinX()] += dist;
                                                upStreamZscores[cr.getBinX()] += zscore;
                                                upStreamVals[cr.getBinX()] += Math.log(1+cr.getCounts());
                                                upStreamCounts[cr.getBinX()]++;

                                                downStreamSums[cr.getBinY()] += (dist * zscore);
                                                downStreamDists[cr.getBinY()] += dist;
                                                downStreamZscores[cr.getBinY()] += zscore;
                                                downStreamVals[cr.getBinY()] += Math.log(1+cr.getCounts());
                                                downStreamCounts[cr.getBinY()]++;
                                                */

                                                upStreamSums[cr.getBinX()] += (dist * (cr.getCounts() + 1));
                                                upStreamVals[cr.getBinX()] += (dist * (poly.getExpectedFromUncompressedBin(dist) + 1));

                                                downStreamSums[cr.getBinY()] += (dist * (cr.getCounts() + 1));
                                                downStreamVals[cr.getBinY()] += (dist * (poly.getExpectedFromUncompressedBin(dist) + 1));

                                            }
                                        }
                                    }
                                }
                            }

                            divide(upStreamSums, upStreamVals);
                            divide(downStreamSums, downStreamVals);

                            /*
                            divide(upStreamSums, upStreamCounts);
                            divide(downStreamSums, downStreamCounts);

                            divide(upStreamZscores, upStreamCounts);
                            divide(downStreamZscores, downStreamCounts);

                            divide(upStreamVals, upStreamCounts);
                            divide(downStreamVals, downStreamCounts);

                            divide(upStreamDists, upStreamCounts);
                            divide(downStreamDists, downStreamCounts);
                            */


                            synchronized (allUpStreamRowSums) {
                                allUpStreamRowSums.put(chrom, upStreamSums);
                            }
                            synchronized (allDownStreamRowSums) {
                                allDownStreamRowSums.put(chrom, downStreamSums);
                            }
                            /*
                            synchronized (allUpStreamZscores) {
                                allUpStreamZscores.put(chrom, upStreamZscores);
                            }
                            synchronized (allDownStreamZscores) {
                                allDownStreamZscores.put(chrom, downStreamZscores);
                            }
                            synchronized (allUpStreamDists) {
                                allUpStreamDists.put(chrom, upStreamDists);
                            }
                            synchronized (allDownStreamDists) {
                                allDownStreamDists.put(chrom, downStreamDists);
                            }
                            synchronized (allUpStreamVals) {
                                allUpStreamVals.put(chrom, upStreamVals);
                            }
                            synchronized (allDownStreamVals) {
                                allDownStreamVals.put(chrom, downStreamVals);
                            }
                            synchronized (allUpStreamCounts) {
                                allUpStreamCounts.put(chrom, upStreamCounts);
                            }
                            synchronized (allDownStreamCounts) {
                                allDownStreamCounts.put(chrom, downStreamCounts);
                            }

 */

                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                        }
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });

        System.out.println("Exporting anchor results...");
        try {
            SeerUtils.exportRowFloatsToBedgraph(allUpStreamRowSums, outputPath + ".upstream.sums.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownStreamRowSums, outputPath + ".downstream.sums.bedgraph", resolution);

            /*
            SeerUtils.exportRowFloatsToBedgraph(allUpStreamZscores, outputPath + ".upstream.zscores.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownStreamZscores, outputPath + ".downstream.zscores.bedgraph", resolution);
            SeerUtils.exportRowIntsToBedgraph(allUpStreamDists, outputPath + ".upstream.dists.bedgraph", resolution);
            SeerUtils.exportRowIntsToBedgraph(allDownStreamDists, outputPath + ".downstream.dists.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allUpStreamVals, outputPath + ".upstream.values.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownStreamVals, outputPath + ".downstream.values.bedgraph", resolution);
            SeerUtils.exportRowIntsToBedgraph(allUpStreamCounts, outputPath + ".upstream.counts.bedgraph", resolution);
            SeerUtils.exportRowIntsToBedgraph(allDownStreamCounts, outputPath + ".downstream.counts.bedgraph", resolution);

             */

        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Anchor strengths complete");
    }

    private Chromosome[] getChromosomes(ChromosomeHandler handler) {
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        if (chrom != null && chrom.length() > 0) {
            String[] chroms = chrom.split(",");
            chromosomes = new Chromosome[chroms.length];
            for (int i = 0; i < chroms.length; i++) {
                chromosomes[i] = handler.getChromosomeFromName(chroms[i]);
            }
        }
        return chromosomes;
    }

    private void divide(float[] sums, int[] counts) {
        for (int k = 0; k < sums.length; k++) {
            if (counts[k] > 0) {
                sums[k] /= counts[k];
            }
        }
    }

    private void divide(float[] sums, float[] counts) {
        for (int k = 0; k < sums.length; k++) {
            if (counts[k] > 0) {
                sums[k] /= counts[k];
            }
        }
    }

    private void divide(int[] sums, int[] counts) {
        for (int k = 0; k < sums.length; k++) {
            if (counts[k] > 0) {
                sums[k] /= counts[k];
            }
        }
    }
}