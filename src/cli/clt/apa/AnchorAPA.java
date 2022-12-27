package cli.clt.apa;

import cli.clt.CommandLineParser;
import cli.utils.apa.APAUtils;
import cli.utils.apa.AnchorAPAScore;
import cli.utils.apa.MultiAPAManager;
import cli.utils.data.SparseContactMatrixWithMasking;
import cli.utils.flags.RegionConfiguration;
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
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorAPA {
    public static String usage = "anchor-apa [--ag-norm] [-k NORM] [--window val]" +
            " [--min-dist val] [--max-dist val] [-r resolution]" +
            " <input.hic> <loops.bedpe> <outfolder>";
    protected final NormalizationType norm;
    protected final int window;
    protected final int matrixWidthL;
    protected final int resolution;
    private final String loopListPath;
    private final File outputDirectory;
    private final Dataset ds;
    private final boolean useAgNorm;

    private final AtomicInteger currChromPair = new AtomicInteger(0);
    private final AtomicInteger currNumLoops = new AtomicInteger(0);
    private final NormalizationType vcNorm = NormalizationHandler.VC;
    private final HiCZoom zoom;
    private final ChromosomeHandler handler;
    private final Feature2DList loopList;
    private final int numTotalLoops;

    public AnchorAPA(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(5000);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 50);
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
        int minPeakDist = parser.getMinDistVal(2 * window);
        int maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);

        matrixWidthL = 2 * window + 1;
        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();

        // needs to be last step
        loopList = loadUniqueDistanceFilteredLoops(handler, minPeakDist, maxPeakDist);
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
        System.out.println("Processing Anchor APA for resolution " + resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new HashMap<>();
        int pairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        List<AnchorAPAScore> scores = new ArrayList<>(numTotalLoops);

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < pairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();

                Matrix matrix = ds.getMatrix(chr1, chr1);
                if (matrix != null) {
                    List<Feature2D> loops = loopList.get(chr1.getIndex(), chr1.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {

                                double[] vector = ds.getNormalizationVector(chr1.getIndex(), zoom,
                                        vcNorm).getData().getValues().get(0);

                                SparseContactMatrixWithMasking scm = new SparseContactMatrixWithMasking(zd,
                                        loops, resolution, window, matrixWidthL, norm);

                                MultiAPAManager manager = new MultiAPAManager(loops, window, resolution,
                                        matrixWidthL, scm, vector, useAgNorm);

                                List<AnchorAPAScore> threadScores = manager.getAnchorAPAScores(chr1, resolution);
                                synchronized (scores) {
                                    scores.addAll(threadScores);
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
        });

        System.out.println("Exporting Anchor-APA results...");
        export(scores);
        System.out.println("Anchor-APA complete");
    }

    private void export(List<AnchorAPAScore> scores) {
        File outPath = new File(outputDirectory, "anchor-apa.bed");
        try (PrintWriter writer = new PrintWriter(outPath.getPath())) {
            writer.println("chrom\tx1\tx2\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb");
            for (AnchorAPAScore score : scores) {
                writer.println(score.getLineForBEDFile());
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error exporting Anchor-APA results: " + e.getMessage());
        }

        outPath = new File(outputDirectory, "anchor-apa.bedgraph");
        try (PrintWriter writer = new PrintWriter(outPath.getPath())) {
            writer.println("chrom\tx1\tx2\tscore");
            for (AnchorAPAScore score : scores) {
                writer.println(score.getLineForBedgraphFile());
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error exporting Anchor-APA results: " + e.getMessage());
        }
    }

    private Feature2DList loadUniqueDistanceFilteredLoops(ChromosomeHandler handler, int minPeakDist, int maxPeakDist) {
        return Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> APAUtils.filterFeaturesBySize(new ArrayList<>(new HashSet<>(features)),
                        minPeakDist, maxPeakDist, resolution), false);
    }
}
