package cli.clt.misc;

import cli.clt.CommandLineParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class NormHack {
    public static String usage = "hack [--res int] <out_folder> <file1.hic,file2.hic,...> <name1,name2,...>";

    public static void run(String[] args, CommandLineParser parser) {
        File outFolder = UNIXTools.makeDir(new File(args[1]));
        String[] filepaths = args[2].split(",");
        String[] names = args[3].split(",");

        Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[0], false, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        NormalizationType norm = parser.getNormOrDefaultScale(ds);
        ds.clearCache(false);
        ds = null;

        int numRows = filepaths.length;
        System.out.println("Number of files: " + numRows);

        int resolution = parser.getResolutionOption(25000);
        if (resolution < 1) {
            resolution = 25000;
        }

        System.out.println("Using resolution: " + resolution);

        float[][] result1 = getAllTheNorms(filepaths, names, handler, resolution, norm);
        MatrixTools.saveMatrixTextNumpy((new File(outFolder, "full_norms.npy")).getAbsolutePath(),
                result1);

        float[][] result2 = cleanupNorms(result1);
        result1 = null;
        MatrixTools.saveMatrixTextNumpy((new File(outFolder, "clean_log_norms.npy")).getAbsolutePath(),
                result2);

        System.out.println("norm hack complete");
    }

    private static float[][] getAllTheNorms(String[] filepaths, String[] names, ChromosomeHandler handler, int resolution,
                                            NormalizationType norm) {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        int[] breakdowns = getBreakdowns(chromosomes, resolution);
        int numColumns = breakdowns[breakdowns.length - 1];

        float[][] data = new float[filepaths.length][numColumns];
        for (int p = 0; p < filepaths.length; p++) {
            final Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[p], false, false,
                    true);
            String prefix = names[p];
            for (int c = 0; c < chromosomes.length; c++) {
                populateNormVector(data[p], breakdowns[c], chromosomes[c], resolution, norm, ds, prefix);
            }
            ds.clearCache(false);
        }

        return data;
    }

    private static void populateNormVector(float[] normValues, int offset, Chromosome chromosome,
                                           int resolution, NormalizationType norm, Dataset ds, String prefix) {
        NormalizationVector nv1 = ds.getNormalizationVector(chromosome.getIndex(), new HiCZoom(resolution), norm);
        if (nv1 != null && nv1.getData() != null && nv1.getData().getValues().size() == 1) {
            int width = (int) ((chromosome.getLength() / resolution) + 1);
            float[] row = toFloats(nv1.getData().getValues().get(0));
            System.arraycopy(row, 0, normValues, offset, Math.min(row.length, width));
        } else {
            System.err.println("Error with " + norm.getLabel() + " for " + prefix +
                    " at chromosome: " + chromosome.getName() + " resolution: " + resolution);
        }
    }

    private static float[] toFloats(double[] doubles) {
        float[] floats = new float[doubles.length];
        for (int k = 0; k < floats.length; k++) {
            floats[k] = (float) doubles[k];
        }
        return floats;
    }

    private static int[] getBreakdowns(Chromosome[] chromosomes, int resolution) {
        int[] bounds = new int[chromosomes.length + 1];
        for (int k = 0; k < chromosomes.length; k++) {
            int width = (int) ((chromosomes[k].getLength() / resolution) + 1);
            bounds[k + 1] = bounds[k] + width;
        }
        return bounds;
    }

    private static float[][] cleanupNorms(float[][] initial) {
        float[][] transposed = MatrixTools.transpose(initial);
        List<Integer> whatToSave = getGoodRows(transposed);
        float[][] filtered = new float[whatToSave.size()][initial.length];
        for (int i = 0; i < whatToSave.size(); i++) {
            System.arraycopy(transposed[whatToSave.get(i)], 0, filtered[i], 0, initial.length);
        }
        transposed = null;
        whatToSave.clear();
        inPlaceLog(filtered);
        return MatrixTools.transpose(filtered);
    }

    private static void inPlaceLog(float[][] filtered) {
        for (int r = 0; r < filtered.length; r++) {
            for (int c = 0; c < filtered[r].length; c++) {
                filtered[r][c] = (float) Math.log(filtered[r][c]);
            }
        }
    }

    private static List<Integer> getGoodRows(float[][] matrix) {
        List<Integer> indices = new ArrayList<>();
        for (int r = 0; r < matrix.length; r++) {
            if (isGood(matrix[r])) {
                indices.add(r);
            }
        }
        return indices;
    }

    private static boolean isGood(float[] row) {
        for (float val : row) {
            if (val > 0) {
                // intentionally using > 0 because it does
                // both a nan and positive check
            } else {
                return false;
            }
        }
        return true;
    }
}
