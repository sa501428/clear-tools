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

package cli.clt.stripes;

import cli.clt.CommandLineParser;
import cli.utils.StrawUtils;
import cli.utils.data.SparseFilteredOEMap;
import cli.utils.seer.SeerUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class Slash {


    // FISH - Fast Identification of Stripes in Hi-C
    // SLASH - Statistical Localization and Annotation of Stripes in Hi-C
    public static String usage = "slash [-k NORM] [-c chrom]" +
            "[--min-dist val] [--max-dist val] [-r resolution] <input.hic> <outfile.bedpe>\n" +
            "find stripes in a Hi-C map";
    private final String outputFile;
    private final Dataset ds;
    private final int minPeakDist, maxPeakDist; // distance between two bins, can be changed in opts
    private final int resolution;
    private final String chrom;
    private final int minLengthStripe = 10;
    private NormalizationType norm = NormalizationHandler.VC;

    public Slash(String[] args, CommandLineParser parser, String name) {
        if (args.length != 3) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(2000);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 50);
        outputFile = args[2];
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
        maxPeakDist = parser.getMaxDistVal(8000000) / resolution;
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }


    public void run() {
        System.out.println("Calling stripes for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(resolution);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        Chromosome[] chromosomes = StrawUtils.getChromosomes(handler, chrom);

        final AtomicInteger currChromPair = new AtomicInteger(0);

        final Feature2DList stripes = new Feature2DList();
        final Feature2DList horizontalStripes = new Feature2DList();
        final Feature2DList verticalStripes = new Feature2DList();
        final Map<Chromosome, float[]> horizontals = new HashMap<>();
        final Map<Chromosome, float[]> verticals = new HashMap<>();

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomes.length) {
                Chromosome chrom = chromosomes[threadPair];
                Matrix matrix = ds.getMatrix(chrom, chrom);
                if (matrix != null) {

                    MatrixZoomData zd = matrix.getZoomData(zoom);
                    if (zd != null) {
                        try {

                            System.out.println("creating dataset for " + chrom.getName());
                            SparseFilteredOEMap map = new SparseFilteredOEMap(ds, zd, norm, chrom, resolution,
                                    minPeakDist, maxPeakDist, zoom, minLengthStripe);
                            System.out.println("Getting horizontal stripes for " + chrom.getName());
                            List<Feature2D> horizontalStripesForChrom = map.getHorizontalStripes();
                            System.out.println("Getting vertical stripes for " + chrom.getName());
                            List<Feature2D> verticalStripesForChrom = map.getVerticalStripes();

                            synchronized (horizontals) {
                                horizontals.put(chrom, map.getHorizontalSignal());
                            }
                            synchronized (verticals) {
                                verticals.put(chrom, map.getVerticalSignal());
                            }

                            map.clear();
                            map = null;

                            synchronized (stripes) {
                                stripes.addByKey(Feature2DList.getKey(chrom, chrom), horizontalStripesForChrom);
                                stripes.addByKey(Feature2DList.getKey(chrom, chrom), verticalStripesForChrom);
                            }
                            synchronized (horizontalStripes) {
                                horizontalStripes.addByKey(Feature2DList.getKey(chrom, chrom), horizontalStripesForChrom);
                            }
                            synchronized (verticalStripes) {
                                verticalStripes.addByKey(Feature2DList.getKey(chrom, chrom), verticalStripesForChrom);
                            }
                            horizontalStripesForChrom.clear();
                            verticalStripesForChrom.clear();


                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                        }
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });

        System.out.println("Exporting results...");
        try {
            SeerUtils.exportRowFloatsToBedgraph(horizontals, outputFile + ".up.bedgraph", resolution);
            SeerUtils.exportRowFloatsToBedgraph(verticals, outputFile + ".down.bedgraph", resolution);
        } catch (IOException e) {
            e.printStackTrace();
        }
        stripes.exportFeatureList(new File(outputFile), false, Feature2DList.ListFormat.NA);
        horizontalStripes.exportFeatureList(new File(outputFile.replace(".bedpe", ".horizontal.bedpe")), false, Feature2DList.ListFormat.NA);
        verticalStripes.exportFeatureList(new File(outputFile.replace(".bedpe", ".vertical.bedpe")), false, Feature2DList.ListFormat.NA);
        System.out.println("SLASH complete");
    }
}