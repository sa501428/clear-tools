package cli.clt;

import cli.Main;
import cli.utils.IGVUtils;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
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
    private final int resolution;

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
    }

    public void run() {
        try {
            aggregate(inputBigWig, bedFile, handler, resolution, outFile, genome);
        } catch (IOException | TribbleIndexNotFoundException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    private static void aggregate(String inputBigWig, String bedFile, ChromosomeHandler handler, int resolution, String outFile, String genomeID) throws IOException, TribbleIndexNotFoundException {
        DataTrack igvTrack = IGVUtils.loadBigWig(inputBigWig, handler);

        Chromosome chrom = handler.getChromosomeFromName("chr1");

        List<LocusScore> loci = IGVUtils.getLocusScores(igvTrack, chrom, 0, 100, resolution,
                WindowFunction.count);
        for (LocusScore locus : loci) {
            System.out.println(locus.getStart() + " " + locus.getEnd() + " " + locus.getScore());
        }


        TribbleFeatureSource src = IGVUtils.loadBed(bedFile, handler);

        Iterator<?> iter = src.getFeatures("chr10", 0, (int) handler.getChromosomeFromName("chr10").getLength());
        while (iter.hasNext()) {
            IGVFeature feature = (IGVFeature) iter.next();
            System.out.println(feature.getStart() + " " + feature.getEnd() + " " + feature.getScore());
        }


    }
}
