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
import cli.utils.anchors.AnchorPeakFinder;
import cli.utils.general.ArrayTools;
import cli.utils.general.BedTools;
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

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorStrength {
    public static String usage = "anchor-strength[-sqrt] [-k NORM] [-c chrom]" +
            "[--min-dist val] [--max-dist val] [-r resolution] <input.hic> <output.stem>\n" +
            "calculate localized row sums near loops";
    private final String outputPath;
    private final Dataset ds;
    private final int minPeakDist, maxPeakDist; // distance between two bins, can be changed in opts
    private final int resolution;
    private NormalizationType norm = NormalizationHandler.VC;
    private final String chrom;
    private final boolean useSqrtScaling;
    private final boolean useExperimentalVersion;

    public AnchorStrength(String[] args, CommandLineParser parser, String name) {
        if (args.length != 3) {
            printUsageAndExit();
        }

        useSqrtScaling = name.contains("sqrt");
        useExperimentalVersion = name.contains("dev");

        resolution = parser.getResolutionOption(2000);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 50);
        outputPath = args[2];
        chrom = parser.getChromosomeOption();

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            if (possibleNorm != null && possibleNorm.length() > 0) {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        System.out.println("Using normalization: " + norm.getLabel());
        minPeakDist = parser.getMinDistVal(20000) / resolution;
        maxPeakDist = parser.getMaxDistVal(5000000) / resolution;
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

        final Map<Chromosome, float[]> allUpStreamOEProd = new HashMap<>();
        final Map<Chromosome, float[]> allDownStreamOEProd = new HashMap<>();
        final Map<Chromosome, float[]> allBothStreamOEProd = new HashMap<>();

        /*
        final Map<Chromosome, float[]> allUpSumProbs = new HashMap<>();
        final Map<Chromosome, float[]> allDownSumProbs = new HashMap<>();
        final Map<Chromosome, float[]> allUpSumProbsTimesDist = new HashMap<>();
        final Map<Chromosome, float[]> allDownSumProbsTimesDist = new HashMap<>();
        final Map<Chromosome, float[]> allUpScaledPercTotal = new HashMap<>();
        final Map<Chromosome, float[]> allDownScaledPercTotal = new HashMap<>();
        */

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
                            float[] upStreamOEP = new float[numEntries];
                            float[] downStreamOEP = new float[numEntries];
                            float[] bothStreamOEP = new float[numEntries];
                            Arrays.fill(upStreamOEP, 1);
                            Arrays.fill(downStreamOEP, 1);
                            Arrays.fill(bothStreamOEP, 1);

                            int[] countsUpstream = new int[numEntries];
                            int[] countsDownstream = new int[numEntries];

                            /*
                            float[] upSumProbs = new float[numEntries];
                            float[] downSumProbs = new float[numEntries];
                            float[] upSumProbsTimesDist = new float[numEntries];
                            float[] downSumProbsTimesDist = new float[numEntries];
                            int[] upDistancesTotal = new int[numEntries];
                            int[] downDistancesTotal = new int[numEntries];
                            */

                            // max p and max d for each locus?
                            // max p x d and max d for each locus?

                            //sum all (d x p) / (sum d)

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
                                            //float zscore = (float) poly.getZscoreForObservedUncompressedBin(dist, cr.getCounts());
                                            if (oe > 2) { // zscore > 1 && oe > 2

                                                upStreamOEP[cr.getBinX()] *= oe;
                                                downStreamOEP[cr.getBinY()] *= oe;

                                                bothStreamOEP[cr.getBinX()] *= oe;
                                                bothStreamOEP[cr.getBinY()] *= oe;

                                                countsUpstream[cr.getBinX()]++;
                                                countsDownstream[cr.getBinY()]++;

                                                /*
                                                double perc = poly.getPercentContact(cr);
                                                if (perc > 0) {
                                                    // total percent up/downstream interaction sum(p)
                                                    upSumProbs[cr.getBinX()] += perc;
                                                    downSumProbs[cr.getBinY()] += perc;

                                                    //distance of influence
                                                    //proportional distance effect sum(d x p)
                                                    upSumProbsTimesDist[cr.getBinX()] += (perc * dist);
                                                    downSumProbsTimesDist[cr.getBinY()] += (perc * dist);

                                                    upDistancesTotal[cr.getBinX()] += dist;
                                                    downDistancesTotal[cr.getBinY()] += dist;
                                                }
                                                */
                                            }
                                        }
                                    }
                                }
                            }

                            //int numLoopyEntries = VectorCleaner.getPercentile(counts, 50, 2);

                            /*
                            float[] upScaledPerc = new float[numEntries];
                            float[] downScaledPerc = new float[numEntries];

                            for (int z = 0; z < upSumProbsTimesDist.length; z++) {
                                if (upDistancesTotal[z] > 0) {
                                    upScaledPerc[z] = upSumProbsTimesDist[z] / upDistancesTotal[z];
                                }
                                if (downDistancesTotal[z] > 0) {
                                    downScaledPerc[z] = downSumProbsTimesDist[z] / downDistancesTotal[z];
                                }

                                upSumProbsTimesDist[z] *= resolution;
                                downSumProbsTimesDist[z] *= resolution;
                            }
                            */

                            for (int z = 0; z < countsUpstream.length; z++) {
                                if (countsUpstream[z] > 0) {
                                    upStreamOEP[z] = (float) (Math.pow(upStreamOEP[z], 1.0 / countsUpstream[z]) * scaleBy(countsUpstream[z]));
                                }
                                if (countsDownstream[z] > 0) {
                                    downStreamOEP[z] = (float) (Math.pow(downStreamOEP[z], 1.0 / countsDownstream[z]) * scaleBy(countsDownstream[z]));
                                }
                                if (countsUpstream[z] + countsDownstream[z] > 0) {
                                    bothStreamOEP[z] = (float) (Math.pow(bothStreamOEP[z], 1.0 / (countsUpstream[z] + countsDownstream[z]))
                                            * scaleBy(countsUpstream[z] + countsDownstream[z]));
                                }
                            }

                            synchronized (allUpStreamOEProd) {
                                allUpStreamOEProd.put(chrom, upStreamOEP);
                            }
                            synchronized (allDownStreamOEProd) {
                                allDownStreamOEProd.put(chrom, downStreamOEP);
                            }
                            synchronized (allBothStreamOEProd) {
                                allBothStreamOEProd.put(chrom, bothStreamOEP);
                            }
                            /*
                            synchronized (allUpSumProbs) {
                                allUpSumProbs.put(chrom, upSumProbs);
                            }
                            synchronized (allDownSumProbs) {
                                allDownSumProbs.put(chrom, downSumProbs);
                            }
                            synchronized (allUpSumProbsTimesDist) {
                                allUpSumProbsTimesDist.put(chrom, upSumProbsTimesDist);
                            }
                            synchronized (allDownSumProbsTimesDist) {
                                allDownSumProbsTimesDist.put(chrom, downSumProbsTimesDist);
                            }
                            synchronized (allUpScaledPercTotal) {
                                allUpScaledPercTotal.put(chrom, upScaledPerc);
                            }
                            synchronized (allDownScaledPercTotal) {
                                allDownScaledPercTotal.put(chrom, downScaledPerc);
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
            SeerUtils.exportRowFloatsToBedgraph(allUpStreamOEProd, outputPath + ".forward.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownStreamOEProd, outputPath + ".reverse.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allBothStreamOEProd, outputPath + ".mixed.bedgraph", resolution);

            /*
            SeerUtils.exportRowFloatsToBedgraph(allUpSumProbs, outputPath + ".upSumProbs.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownSumProbs, outputPath + ".downSumProbs.bedgraph", resolution);

            SeerUtils.exportRowFloatsToBedgraph(allUpSumProbsTimesDist, outputPath + ".upSumProbsTimesDist.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownSumProbsTimesDist, outputPath + ".downSumProbsTimesDist.bedgraph", resolution);

            SeerUtils.exportRowFloatsToBedgraph(allUpScaledPercTotal, outputPath + ".upScaledPercTotal.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(allDownScaledPercTotal, outputPath + ".downScaledPercTotal.bedgraph", resolution);
            */

            BedTools.exportBedFile(new File(outputPath + ".anchors.bed"),
                    AnchorPeakFinder.getPeaks(resolution, allUpStreamOEProd, allDownStreamOEProd, allBothStreamOEProd));


        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Anchor strengths complete");
    }

    private double scaleBy(int n) {
        if (useSqrtScaling) {
            return Math.sqrt(n);
        } else {
            return n;
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
}