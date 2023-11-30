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

package cli.clt.anchor;

import cli.clt.CommandLineParser;
import cli.utils.general.ArrayTools;
import cli.utils.general.MapNMS;
import cli.utils.general.QuadContactRecord;
import cli.utils.general.VectorCleaner;
import cli.utils.general.function.OELogMedianFunction;
import cli.utils.general.function.ScaledGeoMeanFunction;
import cli.utils.general.function.SimpleDivisionFunction;
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

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorStrength {
    private static final float Z_TOP_TEN = 1.28f;
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
    private final int compressionScalar = 5;


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

        final Map<Chromosome, float[][]> allUpMatrices = new HashMap<>();
        final Map<Chromosome, float[][]> allDownMatrices = new HashMap<>();

        final Map<Chromosome, float[][]> allUpMatricesWithNMS = new HashMap<>();
        final Map<Chromosome, float[][]> allDownMatricesWithNMS = new HashMap<>();

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

                            float[][] upMatrix = new float[MapNMS.NUM_ROWS][numEntries];
                            float[][] downMatrix = new float[MapNMS.NUM_ROWS][numEntries];
                            MapNMS.initOERows(upMatrix, downMatrix);

                            LogExpectedZscoreSpline poly = new LogExpectedZscoreSpline(zd, norm, chrom, resolution);

                            double[] vector = ArrayTools.copy(ds.getNormalizationVector(chrom.getIndex(), zoom,
                                    NormalizationHandler.VC).getData().getValues().get(0));
                            VectorCleaner.inPlaceZscore(vector);

                            Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap = new HashMap<>();

                            Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
                            while (it.hasNext()) {
                                ContactRecord cr = it.next();
                                if (cr.getCounts() > 0) {
                                    if (vector[cr.getBinX()] > -1 && vector[cr.getBinY()] > -1) {

                                        int dist = ExpectedUtils.getDist(cr);
                                        if (dist > minPeakDist && dist < maxPeakDist) {

                                            float oe = (float) ((cr.getCounts() + 1) / (poly.getExpectedFromUncompressedBin(dist) + 1));
                                            float zscore = (float) poly.getZscoreForObservedUncompressedBin(dist, cr.getCounts());
                                            float perc = poly.getPercentContact(cr);
                                            if (perc > 0 && (oe > 2 || zscore > Z_TOP_TEN)) { // zscore > 1 && oe > 2

                                                MapNMS.fillInScoreRows(upMatrix, downMatrix,
                                                        cr.getBinX(), cr.getBinY(),
                                                        cr.getCounts(), oe, perc, zscore);

                                                populateMap(contactMap, cr, compressionScalar, oe, perc, zscore);
                                            }
                                        }
                                    }
                                }
                            }

                            float[][] upMatrixWithNMS = new float[MapNMS.NUM_ROWS][numEntries];
                            float[][] downMatrixWithNMS = new float[MapNMS.NUM_ROWS][numEntries];

                            MapNMS.populateAfterNonMaxSuppression(upMatrixWithNMS, downMatrixWithNMS, contactMap, numEntries);
                            clearAll(contactMap);

                            synchronized (allUpMatrices) {
                                allUpMatrices.put(chrom, upMatrix);
                                allDownMatrices.put(chrom, downMatrix);
                                allUpMatricesWithNMS.put(chrom, upMatrixWithNMS);
                                allDownMatricesWithNMS.put(chrom, downMatrixWithNMS);
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

            export(allUpMatrices, "forward");
            export(allDownMatrices, "reverse");
            export(allUpMatricesWithNMS, "forward.nms");
            export(allDownMatricesWithNMS, "reverse.nms");

            //BedTools.exportBedFile(new File(outputPath + ".anchors.bed"),
            //        AnchorPeakFinder.getPeaks(resolution, countScaledAllUpStreamOEProd, countScaledAllDownStreamOEProd));

        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Anchor strengths complete");
    }

    private void export(Map<Chromosome, float[][]> allMatrices, String stem) throws IOException {

        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".raw." + stem + ".bedgraph", resolution,
                MapNMS.COUNT_INDEX, null);
        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".oe." + stem + ".bedgraph", resolution,
                MapNMS.OE_INDEX, null);
        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".perc." + stem + ".bedgraph", resolution,
                MapNMS.PERC_INDEX, null);

        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".oeLogMed." + stem + ".bedgraph", resolution,
                MapNMS.OE_INDEX, new OELogMedianFunction());
        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".zscore." + stem + ".bedgraph", resolution,
                MapNMS.ZSCORE_INDEX, new SimpleDivisionFunction());
        SeerUtils.exportRowFloatsToBedgraph(allMatrices, outputPath + ".scaleGeomOE." + stem + ".bedgraph", resolution,
                MapNMS.OE_INDEX, new ScaledGeoMeanFunction());
    }

    private void clearAll(Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap) {
        for (Map<Integer, List<QuadContactRecord>> map : contactMap.values()) {
            for (List<QuadContactRecord> list : map.values()) {
                list.clear();
            }
            map.clear();
        }
        contactMap.clear();
    }


    private void populateMap(Map<Integer, Map<Integer, List<QuadContactRecord>>> contactMap, ContactRecord cr,
                             int compressionScalar, float oe, float perc, float zscore) {
        int compressedX = cr.getBinX() / compressionScalar;
        int compressedY = cr.getBinY() / compressionScalar;
        if (!contactMap.containsKey(compressedX)) {
            contactMap.put(compressedX, new HashMap<>());
        }
        if (!contactMap.get(compressedX).containsKey(compressedY)) {
            contactMap.get(compressedX).put(compressedY, new LinkedList<>());
        }
        contactMap.get(compressedX).get(compressedY).add(new QuadContactRecord(cr, oe, perc, zscore));
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