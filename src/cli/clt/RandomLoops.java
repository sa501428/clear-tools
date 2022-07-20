package cli.clt;

import cli.utils.sift.SimpleLocation;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class RandomLoops {

    public static void run(String[] args, CommandLineParser parser) {
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1],
                false, false, true);

        int resolution = parser.getResolutionOption(2500000);

        Map<Chromosome, Set<SimpleLocation>> data = new HashMap<>();
        for (Chromosome chromosome : ds.getChromosomeHandler().getAutosomalChromosomesArray()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome, resolution);
            if (matrix == null) {
                continue;
            }
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) {
                continue;
            }

            Set<SimpleLocation> locations = new HashSet<>();

        }
    }
}
