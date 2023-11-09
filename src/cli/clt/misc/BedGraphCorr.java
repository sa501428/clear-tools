package cli.clt.misc;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.BedGraphParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.MatrixTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class BedGraphCorr {

    public static String usage = "bedgraph-corr[-log][-no-zero] [-r resolution] <genomeID> <output.npy> " +
            "<file1.bedgraph> <file2.bedgraph> ... <fileN.bedgraph>";


    public static void run(String[] args, CommandLineParser parser, String name) {
        boolean doLog = name.contains("log");
        boolean noZero = name.contains("no-zero");

        if (args.length < 5) {
            Main.printGeneralUsageAndExit(67, usage);
        }

        int resolution = parser.getResolutionOption(1000);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        String outputName = args[2];
        List<Map<String, float[]>> files = new ArrayList<>();
        for (int k = 3; k < args.length; k++) {
            files.add(BedGraphParser.parse(handler, args[k], resolution));
        }

        if (doLog) {
            applyLog(files);
        }

        if (noZero) {
            removeZeroes(files);
        }

        process(files, outputName);
    }

    private static void removeZeroes(List<Map<String, float[]>> files) {
        for (Map<String, float[]> file : files) {
            for (float[] values : file.values()) {
                for (int i = 0; i < values.length; i++) {
                    if (values[i] == 0) {
                        values[i] = Float.NaN;
                    }
                }
            }
        }
    }

    private static void applyLog(List<Map<String, float[]>> files) {
        for (Map<String, float[]> file : files) {
            for (float[] values : file.values()) {
                for (int i = 0; i < values.length; i++) {
                    values[i] = (float) Math.log(values[i]);
                }
            }
        }
    }

    private static void process(List<Map<String, float[]>> files, String outputName) {

        if (files.size() == 2) {
            double accuracy = getCorrelation(files.get(0), files.get(1));
            System.out.println("Correlation " + accuracy);
        } else {
            double[][] result = new double[files.size()][files.size()];
            for (int i = 0; i < result.length; i++) {
                result[i][i] = 1;
                for (int j = i + 1; j < result.length; j++) {
                    result[i][j] = getCorrelation(files.get(i), files.get(j));
                    result[j][i] = result[i][j];
                }
            }
            MatrixTools.saveMatrixTextNumpy(outputName, result);
        }
    }

    private static double getCorrelation(Map<String, float[]> file1, Map<String, float[]> file2) {
        double sumX = 0.0;
        double sumY = 0.0;
        int count = 0;

        for (String key : file1.keySet()) {
            if (file2.containsKey(key)) {
                float[] x = file1.get(key);
                float[] y = file2.get(key);
                if (x == null || y == null) {
                    continue;
                }

                for (int i = 0; i < x.length; ++i) {
                    if (isValidComparison(x[i], y[i])) {
                        sumX += x[i];
                        sumY += y[i];
                        count++;
                    }
                }
            }
        }

        double muX = sumX / (double) count;
        double muY = sumY / (double) count;
        double dotProduct = 0.0;
        double normX = 0.0;
        double normY = 0.0;

        for (String key : file1.keySet()) {
            if (file2.containsKey(key)) {
                float[] x = file1.get(key);
                float[] y = file2.get(key);
                if (x == null || y == null) {
                    continue;
                }

                for (int i = 0; i < x.length; ++i) {
                    if (isValidComparison(x[i], y[i])) {
                        double nX = x[i] - muX;
                        double nY = y[i] - muY;
                        dotProduct += nX * nY;
                        normX += nX * nX;
                        normY += nY * nY;
                    }
                }
            }
        }

        return dotProduct / Math.sqrt(normX * normY);
    }

    private static boolean isValidComparison(float x, float v) {
        return (isNotZero(x) || isNotZero(v)) && !Float.isNaN(x) && !Float.isNaN(v);
    }

    private static boolean isNotZero(float x) {
        return !(Math.abs(x) < 1e-10);
    }
}
