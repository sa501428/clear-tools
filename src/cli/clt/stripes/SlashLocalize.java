package cli.clt.stripes;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.data.SparseContactMatrixOfSpecificRegionsOnly;
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
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class SlashLocalize {
    private static final NormalizationType NONE = NormalizationHandler.NONE;
    public static String usage = "slash-localize [--res int] <input.hic> <stripes.bedpe> <outfile.bedpe>";

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        int resolution = parser.getResolutionOption(10);
        int window = parser.getWindowSizeOption(5000 / resolution);

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, false,
                resolution > 50);
        String stripeListPath = args[2];
        String outFile = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();

        Feature2DList stripeList = Feature2DParser.loadFeatures(stripeListPath, handler,
                true, null, false);

        if (Main.printVerboseComments) System.out.println("Number of stripes: " + stripeList.getNumTotalFeatures());

        final Feature2DList horizontalLocalized = new Feature2DList();
        final Feature2DList verticalLocalized = new Feature2DList();
        final Feature2DList allLocalized = new Feature2DList();

        localize(dataset, stripeList, handler, resolution, horizontalLocalized, verticalLocalized, allLocalized, window);
        allLocalized.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
        horizontalLocalized.exportFeatureList(new File(outFile.replaceAll(".bedpe", ".horizontal.bedpe")), false, Feature2DList.ListFormat.NA);
        verticalLocalized.exportFeatureList(new File(outFile.replaceAll(".bedpe", ".vertical.bedpe")), false, Feature2DList.ListFormat.NA);
        if (Main.printVerboseComments) System.out.println("stripe localization complete");
    }

    private static void localize(Dataset ds, Feature2DList stripeList, ChromosomeHandler handler,
                                 int resolution, Feature2DList horizontalLocalized, Feature2DList verticalLocalized,
                                 Feature2DList allLocalized, int compressionSize) {

        HiCZoom zoom = new HiCZoom(resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), true);

        int numTotalLoops = stripeList.getNumTotalFeatures();
        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger currNumDone = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                Matrix matrix = ds.getMatrix(chr1, chr2);
                if (matrix != null) {
                    if (Main.printVerboseComments)
                        System.out.println("Processing " + chr1.getName() + " " + chr2.getName());

                    List<Feature2D> stripes = filterForMinSize(stripeList.get(chr1.getIndex(), chr2.getIndex()),
                            resolution, compressionSize);
                    if (stripes.size() > 0) {
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {
                                if (Main.printVerboseComments)
                                    System.out.println("Num loops to process: " + stripes.size());
                                List<Feature2D> horizontalStripes = new LinkedList<>();
                                List<Feature2D> verticalStripes = new LinkedList<>();

                                SparseContactMatrixOfSpecificRegionsOnly scm = new SparseContactMatrixOfSpecificRegionsOnly(zd,
                                        stripes, resolution, resolution, NONE);

                                if (Main.printVerboseComments) System.out.println("Sparse data loaded");
                                try {

                                    for (Feature2D stripe : stripes) {
                                        int binX1 = (int) (stripe.getStart1() / resolution) - 1;
                                        int binX2 = (int) (stripe.getEnd1() / resolution) + 1;
                                        int binY1 = (int) (stripe.getStart2() / resolution) - 1;
                                        int binY2 = (int) (stripe.getEnd2() / resolution) + 1;

                                        float[][] data = scm.getRegion(binX1, binY1, binX2, binY2);
                                        addToAll(data, getAverage(data));

                                        if (stripe.getWidth2() > stripe.getWidth1()) { // horizontal loop
                                            data = collapseColumns(data, compressionSize);
                                        } else if (stripe.getWidth1() > stripe.getWidth2()) { // vertical loop
                                            data = collapseRows(data, compressionSize);
                                        } else { // should not be happening
                                            data = null;
                                        }

                                        data = normalize(data);

                                        if (stripe.getWidth2() > stripe.getWidth1()) { // horizontal loop
                                            float[] rowSums = getRowSums(data);
                                            int rowIndex = getBestIndex(rowSums);
                                            if (rowIndex > 0 && rowIndex < data.length - 1) {
                                                horizontalStripes.add(makeLocalizedHorizontalStripe(stripe, binX1 + rowIndex, resolution));
                                            }
                                        } else if (stripe.getWidth1() > stripe.getWidth2()) { // vertical loop
                                            float[] colSums = getColSums(data);
                                            int colIndex = getBestIndex(colSums);
                                            if (colIndex > 0 && colIndex < data[0].length - 1) {
                                                verticalStripes.add(makeLocalizedVerticalStripe(stripe, binY1 + colIndex, resolution));
                                            }
                                        }
                                    }
                                } catch (Exception e) {
                                    System.err.println(e.getMessage());
                                    e.printStackTrace();
                                }

                                synchronized (horizontalLocalized) {
                                    horizontalLocalized.addByKey(Feature2DList.getKey(chr1, chr2), horizontalStripes);
                                }
                                synchronized (verticalLocalized) {
                                    verticalLocalized.addByKey(Feature2DList.getKey(chr1, chr2), verticalStripes);
                                }
                                synchronized (allLocalized) {
                                    allLocalized.addByKey(Feature2DList.getKey(chr1, chr2), horizontalStripes);
                                    allLocalized.addByKey(Feature2DList.getKey(chr1, chr2), verticalStripes);
                                }

                                System.out.println(((int) Math.floor((100.0 * currNumDone.get()) / numTotalLoops)) + "% ");

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
    }

    private static int getBestIndex(float[] sums) {

        int bestIndex = -1;
        float bestVal = 0;
        for (int i = 1; i < sums.length - 1; i++) {
            float sum = getSumOfThree(sums, i);
            if (sum > bestVal) {
                bestVal = sum;
                bestIndex = i;
            }
        }

        return bestIndex;
    }

    private static float getSumOfThree(float[] data, int index) {
        return data[index - 1] + data[index] + data[index + 1];
    }

    private static Feature2D makeLocalizedHorizontalStripe(Feature2D stripe, int rowBin, int resolution) {
        long gx = (long) rowBin * resolution;
        return new Feature2D(stripe.getFeatureType(), stripe.getChr1(), gx, gx + resolution,
                stripe.getChr2(), stripe.getStart2(), stripe.getEnd2(), stripe.getColor(), stripe.getAttributes());
    }

    private static Feature2D makeLocalizedVerticalStripe(Feature2D stripe, int colBin, int resolution) {
        long gy = (long) colBin * resolution;
        return new Feature2D(stripe.getFeatureType(), stripe.getChr1(), stripe.getStart1(), stripe.getEnd1(),
                stripe.getChr2(), gy, gy + resolution, stripe.getColor(), stripe.getAttributes());
    }

    private static float[][] normalize(float[][] data) {
        float[] rowSums = getRowSums(data);
        float[] colSums = getColSums(data);

        float[][] normalized = new float[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                normalized[i][j] = (float) (data[i][j] / Math.sqrt(rowSums[i] * colSums[j]));
            }
        }
        return normalized;
    }

    private static float[] getColSums(float[][] data) {
        float[] colSums = new float[data[0].length];
        for (int i = 0; i < data[0].length; i++) {
            for (float[] loop : data) {
                colSums[i] += loop[i];
            }
        }
        return colSums;
    }

    private static float[] getRowSums(float[][] data) {
        float[] rowSums = new float[data.length];
        for (int i = 0; i < data.length; i++) {
            for (float loop : data[i]) {
                rowSums[i] += loop;
            }
        }
        return rowSums;
    }

    private static float[][] collapseRows(float[][] data, int compressionSize) {
        float[][] collapsed = new float[(data.length / compressionSize) + 1][data[0].length];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                collapsed[i / compressionSize][j] += data[i][j];
            }
        }
        return collapsed;
    }

    private static float[][] collapseColumns(float[][] data, int compressionSize) {
        float[][] collapsed = new float[data.length][(data[0].length / compressionSize) + 1];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                collapsed[i][j / compressionSize] += data[i][j];
            }
        }
        return collapsed;
    }

    private static void addToAll(float[][] data, float averageVal) {
        for (float[] loop : data) {
            for (int i = 0; i < loop.length; i++) {
                loop[i] += averageVal;
            }
        }
    }

    private static float getAverage(float[][] data) {
        double sum = 0;
        for (float[] loop : data) {
            for (float loop2 : loop) {
                sum += loop2;
            }
        }
        return (float) (sum / (data.length * data[0].length));
    }

    private static List<Feature2D> filterForMinSize(List<Feature2D> feature2DS, int resolution, int compressionSize) {
        List<Feature2D> filteredBySize = new LinkedList<>();
        if (feature2DS != null) {
            for (Feature2D loop : feature2DS) {
                if (loop.getWidth1() != loop.getWidth2()) {
                    long width = Math.max(loop.getWidth1(), loop.getWidth2());
                    if (compressedSize(width, resolution, compressionSize) > 1) {
                        filteredBySize.add(loop);
                    }
                }
            }
        }
        return filteredBySize;
    }

    private static int compressedSize(long width, int resolution, int compressionSize) {
        return (int) ((width / resolution) / compressionSize);
    }
}
