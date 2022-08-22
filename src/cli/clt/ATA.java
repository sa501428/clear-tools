package cli.clt;

import cli.Main;
import cli.utils.general.IGVTools;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class ATA {

    private final String inputBigWig;
    private final String bedFile;
    private final String outFile;
    private final String genome;
    private final ChromosomeHandler handler;
    private final int resolution, window;

    public static void main(String[] args) {
        // test
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes("hg38");
        String inputBigWig = "/Users/muhammad/Desktop/ATA/S3_w73_lvent_signal.bw";
        DataTrack bigWig = IGVTools.loadBigWig(inputBigWig, handler);

        Chromosome chrom = handler.getChromosomeFromName("chr1");
        int gStart = 1000000;
        int gEnd = 1001000;
        int resolution = 1;


        List<LocusScore> loci = IGVTools.getLocusScores(bigWig, chrom, gStart, gEnd, resolution, WindowFunction.count);
        if (loci.size() > 0) {
            for (LocusScore locus : loci) {
                System.out.println(locus.getStart() + " " + locus.getEnd() + " " + locus.getScore());
            }
        }
    }

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
        final int width = 2 * window + 1;
        final double[] totalAccumulation = new double[width];
        final long[] totalNumber = new long[1];
        final int gWindow = window * resolution;
        final int gWidth = width * resolution;

        TribbleFeatureSource bedFile = IGVTools.loadBed(inputBedFile, handler);
        DataTrack bigWig = IGVTools.loadBigWig(inputBigWig, handler);

        AtomicInteger cIndex = new AtomicInteger(0);

        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        ParallelizationTools.launchParallelizedCode(() -> {
            double[] accumulation = new double[width];
            int peakCounter = 0;

            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {

                try {
                    Chromosome chrom = chromosomes[currIndex];
                    System.out.println("Handling " + chrom.getName());
                    Iterator<?> iter = bedFile.getFeatures(chrom.getName(), 0, (int) chrom.getLength());
                    while (iter.hasNext()) {
                        IGVFeature feature = (IGVFeature) iter.next();
                        //System.out.println(feature.getStart() + " " + feature.getEnd() + " " + feature.getScore());

                        int midPoint = getIntervalMidpoint(feature);
                        int gStart = midPoint - gWindow;
                        int gEnd = gStart + gWidth + 1;

                        List<LocusScore> loci = IGVTools.getLocusScores(bigWig, chrom, gStart, gEnd, 1, WindowFunction.count);
                        if (loci.size() > 0) {
                            peakCounter++;
                            for (LocusScore locus : loci) {
                                if (locus.getScore() > 0) {
                                    int relativeX = (locus.getStart() - gStart) / resolution;
                                    if (relativeX >= 0 && relativeX < width) {
                                        accumulation[relativeX] += locus.getScore();
                                    }
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }

                currIndex = cIndex.getAndIncrement();
            }

            synchronized (totalAccumulation) {
                for (int k = 0; k < totalAccumulation.length; k++) {
                    if (accumulation[k] > 0) {
                        totalAccumulation[k] += accumulation[k];
                    }
                }
                if (peakCounter > 0) {
                    totalNumber[0] += peakCounter;
                }
            }
        });

        normalizeInPlace(totalAccumulation, totalNumber[0]);
        MatrixTools.saveMatrixTextNumpy(outFile + ".npy", totalAccumulation);
    }

    public void run() {
        try {
            aggregate(inputBigWig, bedFile, handler, resolution, outFile, window);
        } catch (IOException | TribbleIndexNotFoundException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    private static void normalizeInPlace(double[] vector, long counts) {
        for (int k = 0; k < vector.length; k++) {
            vector[k] /= counts;
        }
    }

    private static int getIntervalMidpoint(IGVFeature feature) {
        return (feature.getStart() + feature.getEnd()) / 2;
    }
}
