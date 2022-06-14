package cli.clt;

import cli.Main;
import cli.utils.IGVUtils;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.MatrixTools;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

public class ATA {

    private final String inputBigWig;
    private final String bedFile;
    private final String outFile;
    private final String genome;
    private final ChromosomeHandler handler;
    private final int resolution, window;

    public ATA(String[] args, CommandLineParser parser) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        // signal.bw peaks.bed outfile.npy
        inputBigWig = args[1];
        bedFile = args[2];
        outFile = args[3];
        genome = args[4];
        handler = ChromosomeTools.loadChromosomes(genome);
        resolution = parser.getResolutionOption(1);
        window = parser.getWindowSizeOption(100);
    }

    private static void aggregate(String inputBigWig, String inputBedFile, ChromosomeHandler handler, int resolution,
                                  String outFile, int window) throws IOException, TribbleIndexNotFoundException {
        int width = 2 * window + 1;
        double[] accumulation = new double[width];

        TribbleFeatureSource bedFile = IGVUtils.loadBed(inputBedFile, handler);
        DataTrack bigWig = IGVUtils.loadBigWig(inputBigWig, handler);
        int peakCounter = 0;

        for (Chromosome chrom : handler.getAutosomalChromosomesArray()) {
            Iterator<?> iter = bedFile.getFeatures("chr10", 0, (int) handler.getChromosomeFromName("chr10").getLength());
            while (iter.hasNext()) {
                IGVFeature feature = (IGVFeature) iter.next();
                //System.out.println(feature.getStart() + " " + feature.getEnd() + " " + feature.getScore());

                int midPoint = getIntervalMidpoint(feature);
                int gStart = midPoint - window;
                int gEnd = gStart + width + 1;

                List<LocusScore> loci = IGVUtils.getLocusScores(bigWig, chrom, gStart, gEnd, resolution, WindowFunction.count);
                if (loci.size() > 0) {
                    peakCounter++;
                    for (LocusScore locus : loci) {
                        int relativeX = locus.getStart() - gStart;
                        if (relativeX >= 0 && relativeX < width) {
                            accumulation[relativeX] += locus.getScore();
                        }
                    }
                }
            }
        }

        normalizeInPlace(accumulation, peakCounter);
        MatrixTools.saveMatrixTextNumpy(outFile + ".npy", accumulation);
    }

    private static void normalizeInPlace(double[] vector, int peakCounter) {
        for (int k = 0; k < vector.length; k++) {
            vector[k] /= peakCounter;
        }
    }

    private static int getIntervalMidpoint(IGVFeature feature) {
        return (feature.getStart() + feature.getEnd()) / 2;
    }

    public void run() {
        try {
            aggregate(inputBigWig, bedFile, handler, resolution, outFile, window);
        } catch (IOException | TribbleIndexNotFoundException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
}
