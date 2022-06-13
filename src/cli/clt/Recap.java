package cli.clt;

import cli.Main;
import cli.utils.HiCUtils;
import cli.utils.RecapTools;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Recap {

    private NormalizationType norm;

    public Recap(String[] args, CommandLineParser parser) {


        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        // recap <loops.bedpe> <outfolder> <file1.hic,file2.hic,...> <name1,name2,...>
        String loopListPath = args[1];
        File outFolder = UNIXTools.makeDir(new File(args[2]));

        String[] filepaths = args[3].split(",");
        String[] names = args[4].split(",");

        Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[0], false, false,
                true);

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            if (possibleNorm != null && possibleNorm.length() > 0) {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } else {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
            }
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "VC_SQRT", "VC"});
        }
        System.out.println("Using normalization: " + norm.getLabel());

        ds.clearCache(false);
        ds = null;

        int numRows = loopList.getNumTotalFeatures();
        System.out.println("Number of loops: " + numRows);

        int resolution = parser.getResolutionOption(1000);
        if (resolution < 1) {
            resolution = 1000;
        }

        int window = parser.getWindowSizeOption(0);
        if (window < 1) {
            window = 5;
        }

        System.out.println("Using resolution: " + resolution);

        Feature2DList refinedLoops = localize(filepaths, names, loopList, handler, resolution, window, norm);
        refinedLoops.exportFeatureList(new File(outFolder, "recap.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("pinpoint complete");
    }

    private static Feature2DList localize(String[] filepaths, String[] names, Feature2DList loopList,
                                          ChromosomeHandler handler, int resolution, int window,
                                          NormalizationType norm) {

        if (Main.printVerboseComments) {
            System.out.println("Start Recap/Compile process");
        }

        HiCZoom zoom = new HiCZoom(resolution);
        int matrixWidth = 1 + 2 * window;

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final Object key = new Object();

        for (int di = 0; di < filepaths.length; di++) {

            Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[di], false, false,
                    true);
            String prefix = names[di] + "_";


            final AtomicInteger currChromPair = new AtomicInteger(0);
            final AtomicInteger currNumLoops = new AtomicInteger(0);

            ParallelizationTools.launchParallelizedCode(() -> {

                int threadPair = currChromPair.getAndIncrement();
                while (threadPair < chromosomePairCounter) {
                    RegionConfiguration config = chromosomePairs.get(threadPair);
                    Chromosome chrom1 = config.getChr1();
                    Chromosome chrom2 = config.getChr2();

                    List<Feature2D> loops = loopList.get(chrom1.getIndex(), chrom2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom1, chrom2, zoom);

                        if (zd == null) {
                            System.err.println("ZD is null " + chrom1.getName() + "_" + chrom2.getName());
                            System.exit(9);
                        }

                        ExpectedValueFunction df = ds.getExpectedValuesOrExit(zd.getZoom(), norm, chrom1, true, false);

                        if (df == null) {
                            System.err.println("Expected is null " + chrom1.getName() + "_" + chrom2.getName());
                            System.exit(10);
                        }

                        float pseudoCount = getMedianExpectedAt(9000000 / resolution,
                                window, chrom1.getIndex(), df);
                        double superDiagonal = df.getExpectedValue(chrom1.getIndex(), 1);

                        try {
                            for (Feature2D loop : loops) {
                                float[][] obsMatrix = new float[matrixWidth][matrixWidth];
                                float[][] eMatrix = new float[matrixWidth][matrixWidth];

                                Utils.addLocalizedData(obsMatrix, zd, loop, matrixWidth, resolution, window, norm, key);
                                Utils.fillInExpectedMatrix(eMatrix, loop, matrixWidth, df, chrom1.getIndex(),
                                        resolution, window);

                                // MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_raw.npy")).getAbsolutePath(), output);

                                Map<String, String> attributes = RecapTools.getStats(loop, obsMatrix, eMatrix,
                                        resolution, window, df, superDiagonal, pseudoCount);
                                for (String akey : attributes.keySet()) {
                                    loop.addStringAttribute(prefix + akey, attributes.get(akey));
                                }

                                if (currNumLoops.incrementAndGet() % 100 == 0) {
                                    System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                }
                            }
                            zd.clearCache();

                            System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                        }
                    }
                    threadPair = currChromPair.getAndIncrement();
                }
            });
        }

        return loopList;
    }

    private static float getMedianExpectedAt(int d0, int dx, int chromIndex, ExpectedValueFunction df) {
        List<Double> values = new ArrayList<>();
        for (int dist = d0 - dx; dist < d0 + dx; dist++) {
            double expected = df.getExpectedValue(chromIndex, dist);
            if (expected > 0) {
                values.add(expected);
            }
        }
        float median = QuickMedian.fastMedian(values);
        if (median > 0) return median;
        return 0f;
    }
}
