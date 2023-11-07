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

import java.util.Arrays;
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
        minPeakDist = parser.getMinDistVal(Math.max(10, 1000 / resolution));
        maxPeakDist = parser.getMaxDistVal(5000000 / resolution);
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

        final Map<Chromosome, double[]> allUpStreamOEProd = new HashMap<>();
        final Map<Chromosome, double[]> allDownStreamOEProd = new HashMap<>();
        final Map<Chromosome, double[]> allBothStreamOEProd = new HashMap<>();
        final Map<Chromosome, float[]> allBothStreamZscores = new HashMap<>();

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
                            double[] upStreamSums = new double[numEntries];
                            double[] downStreamSums = new double[numEntries];
                            double[] bothStreamSums = new double[numEntries];
                            Arrays.fill(upStreamSums, 1);
                            Arrays.fill(downStreamSums, 1);
                            Arrays.fill(bothStreamSums, 1);

                            float[] bothStreamZscores = new float[numEntries];

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
                                            if (oe > 2) { // zscore > 1 && oe > 2

                                                upStreamSums[cr.getBinX()] *= oe;
                                                downStreamSums[cr.getBinY()] *= oe;

                                                bothStreamSums[cr.getBinX()] *= oe;
                                                bothStreamSums[cr.getBinY()] *= oe;

                                                bothStreamZscores[cr.getBinX()] += zscore;
                                                bothStreamZscores[cr.getBinY()] += zscore;

                                                upStreamCounts[cr.getBinX()]++;
                                                downStreamCounts[cr.getBinY()]++;
                                            }
                                        }
                                    }
                                }
                            }

                            //normalizeByPercentile(upStreamSums, 50);
                            //normalizeByPercentile(downStreamSums, 50);
                            int numLoopyEntries = VectorCleaner.getPercentile(upStreamCounts, 50, 2);
                            takeNthRoot(upStreamSums, numLoopyEntries);
                            takeNthRoot(downStreamSums, numLoopyEntries);
                            takeNthRoot(bothStreamSums, numLoopyEntries);
                            divide(bothStreamZscores, numLoopyEntries);

                            synchronized (allUpStreamOEProd) {
                                allUpStreamOEProd.put(chrom, upStreamSums);
                            }
                            synchronized (allDownStreamOEProd) {
                                allDownStreamOEProd.put(chrom, downStreamSums);
                            }
                            synchronized (allBothStreamOEProd) {
                                allBothStreamOEProd.put(chrom, bothStreamSums);
                            }
                            synchronized (allBothStreamZscores) {
                                allBothStreamZscores.put(chrom, bothStreamZscores);
                            }


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
            SeerUtils.exportRowDoublesToBedgraph(allUpStreamOEProd, outputPath + ".upstream.sums.bedgraph", resolution);
            SeerUtils.exportRowDoublesToBedgraph(allDownStreamOEProd, outputPath + ".downstream.sums.bedgraph", resolution);
            SeerUtils.exportRowDoublesToBedgraph(allBothStreamOEProd, outputPath + ".bothstream.sums.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allBothStreamZscores, outputPath + ".bothstream.zscores.bedgraph", resolution);

            //BedTools.exportBedFile(new File(outputPath + ".disoriented.anchors.bed"), getPeaks(resolution, allBothStreamOEProd, allBothStreamZscores));


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

    private void divide(float[] vector, int scalar) {
        for (int k = 0; k < vector.length; k++) {
            vector[k] /= scalar;
        }
    }

    private void takeNthRoot(double[] vector, int n) {
        for (int k = 0; k < vector.length; k++) {
            vector[k] = Math.pow(vector[k], 1.0 / n);
        }
    }


    private void normalizeByPercentile(double[] vector, int percentile) {
        double perc = VectorCleaner.getPercentile(vector, percentile, 2);
        for (int k = 0; k < vector.length; k++) {
            vector[k] /= perc;
        }
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