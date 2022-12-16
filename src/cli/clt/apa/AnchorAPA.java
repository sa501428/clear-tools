package cli.clt.apa;

import cli.clt.CommandLineParser;
import cli.utils.apa.APAUtils;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
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

public class AnchorAPA {
    private final static int C0 = 0, C1 = 1, R0 = 2, LOOPID = 3;
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
    private final int minPeakDist; // distance between two bins, can be changed in opts
    private final int maxPeakDist;
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
        minPeakDist = parser.getMinDistVal(2 * window);
        maxPeakDist = parser.getMaxDistVal(Integer.MAX_VALUE);

        matrixWidthL = 2 * window + 1;
        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();

        // needs to be last step
        loopList = loadUniqueDistanceFilteredLoops(handler);
        numTotalLoops = loopList.getNumTotalFeatures();
        if (numTotalLoops < 1) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }
    }

    public static List<float[]> initArrays(int numLoops, int matrixSize) {
        List<float[]> list = new ArrayList<>(numLoops);
        for (int i = 0; i < numLoops; i++) {
            list.add(new float[matrixSize]);
        }
        return list;
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

                                List<float[][]> outputs = initList(loops.size(), matrixWidthL);
                                processLoopsForRegion(zd, loops, outputs, currNumLoops);

                                if (useAgNorm) {
                                    List<float[]> rowSums = initArrays(loops.size(), matrixWidthL);
                                    List<float[]> colSums = initArrays(loops.size(), matrixWidthL);
                                    doAggregateNormalization(chr1, zoom, vcNorm, loops, rowSums, colSums);
                                }

                                // todo get the APA score per anchor
                                // todo get OE per anchor??

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

        System.out.println("Exporting APA results...");
        // todo
        // export();
        System.out.println("APA complete");
    }

    private void doAggregateNormalization(Chromosome chr1, HiCZoom zoom,
                                          NormalizationType vcNorm, List<Feature2D> loops,
                                          List<float[]> rowSums, List<float[]> colSums) {
        double[] vector = ds.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);

        for (int k = 0; k < loops.size(); k++) {
            int binXStart = (int) ((loops.get(k).getMidPt1() / resolution) - window);
            int binYStart = (int) ((loops.get(k).getMidPt2() / resolution) - window);
            APAUtils.addLocalSums(rowSums.get(k), vector, binXStart);
            APAUtils.addLocalSums(colSums.get(k), vector, binYStart);
        }
    }

    private Feature2DList loadUniqueDistanceFilteredLoops(ChromosomeHandler handler) {
        return Feature2DParser.loadFeatures(loopListPath, handler, false,
                (chr, features) -> {
                    return APAUtils.filterFeaturesBySize(new ArrayList<>(new HashSet<>(features)),
                            minPeakDist, maxPeakDist, resolution);
                }, false);
    }

    protected void processLoopsForRegion(MatrixZoomData zd, List<Feature2D> loops,
                                         List<float[][]> outputs, AtomicInteger currNumLoops) {
        Map<Integer, List<int[]>> loopListAsMap = convertToMap(loops);
        Set<Integer> allColumns = getColumns(loops);
        int counter = 0;

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (loopListAsMap.containsKey(cr.getBinX())) {
                    if (allColumns.contains(cr.getBinY())) {
                        populateMatrixIfApplicable(outputs, cr, loopListAsMap.get(cr.getBinX()));
                        if (counter++ % 1000 == 0) {
                            System.out.print(".");
                        }
                    }
                }
            }
        }

        currNumLoops.addAndGet(loops.size());

        loopListAsMap.clear();
        allColumns.clear();
        loopListAsMap = null;
        allColumns = null;
    }

    private void populateMatrixIfApplicable(List<float[][]> matrices, ContactRecord cr, List<int[]> allBounds) {
        for (int[] bounds : allBounds) {
            if (cr.getBinY() >= bounds[C0] && cr.getBinY() < bounds[C1]) {
                int relativeX = cr.getBinX() - bounds[R0];
                int relativeY = cr.getBinY() - bounds[C0];
                matrices.get(bounds[LOOPID])[relativeX][relativeY] += cr.getCounts();
            }
        }
    }

    private Set<Integer> getColumns(List<Feature2D> loops) {
        Set<Integer> columns = new HashSet<>();
        for (Feature2D loop : loops) {
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binYEnd = binYStart + matrixWidthL;

            for (int c = binYStart; c < binYEnd; c++) {
                columns.add(c);
            }
        }
        return columns;
    }

    private Map<Integer, List<int[]>> convertToMap(List<Feature2D> loops) {
        Map<Integer, List<int[]>> loopListAsMap = new HashMap<>();
        for (int k = 0; k < loops.size(); k++) {
            int binXStart = (int) ((loops.get(k).getMidPt1() / resolution) - window);
            int binYStart = (int) ((loops.get(k).getMidPt2() / resolution) - window);
            int binXEnd = binXStart + matrixWidthL;
            int binYEnd = binYStart + matrixWidthL;

            for (int r = binXStart; r < binXEnd; r++) {
                if (!loopListAsMap.containsKey(r)) {
                    loopListAsMap.put(r, new LinkedList<>());
                }
                // C0 = 0, C1 = 1, R0 = 2, LOOPID
                loopListAsMap.get(r).add(new int[]{binYStart, binYEnd, binXStart, k});
            }
        }
        return loopListAsMap;
    }

    public static List<float[][]> initList(int numLoops, int matrixSize) {
        List<float[][]> list = new ArrayList<>(numLoops);
        for (int i = 0; i < numLoops; i++) {
            list.add(new float[matrixSize][matrixSize]);
        }
        return list;
    }
}
