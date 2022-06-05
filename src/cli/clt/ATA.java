package cli.clt;

import cli.Main;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.track.*;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class ATA {

    private final static double log2 = Math.log(2.0D);

    public static void run(String[] args, int resolution, boolean exportNPY) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        // signal.bw peaks.bed outfile.npy

        String inputBigWig = args[1];
        String bedFile = args[2];
        String outFile = args[3];
        String genome = args[4];

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genome);

        aggregate(inputBigWig, bedFile, handler, resolution, outFile, genome);
    }

    public static void main(String[] args) throws IOException, TribbleIndexNotFoundException {
        String inputBed = "https://www.encodeproject.org/files/ENCFF073ORT/@@download/ENCFF073ORT.bed.gz";
        ResourceLocator locator = new ResourceLocator(inputBed);

        String genomeID = "hg38";
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        List<Chromosome> chrList = new ArrayList<>(Arrays.asList(handler.getChromosomeArray()));
        List<org.broad.igv.feature.Chromosome> igvChrList = new ArrayList<>();
        for (Chromosome chrom : chrList) {
            igvChrList.add(chrom.toIGVChromosome());
        }

        Genome genome = new Genome(genomeID, igvChrList);

        TribbleFeatureSource src = TribbleFeatureSource.getFeatureSource(locator, genome);

        Iterator<?> iter = src.getFeatures("chr1", 180000, 182000);
        while (iter.hasNext()) {
            IGVFeature feature = (IGVFeature) iter.next();
            System.out.println(feature.getStart() + " " + feature.getEnd() + " " + feature.getScore());
        }
    }

    public static void main2(String[] args) {

        String inputBigWig = "https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/loops/S5_PC3_signal.bw";
        String genomeID = "hg38";
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        List<Chromosome> chrList = new ArrayList<>(Arrays.asList(handler.getChromosomeArray()));
        List<org.broad.igv.feature.Chromosome> igvChrList = new ArrayList<>();
        for (Chromosome chrom : chrList) {
            igvChrList.add(chrom.toIGVChromosome());
        }

        Genome genome = new Genome(genomeID, igvChrList);
        ResourceLocator locator = new ResourceLocator(inputBigWig);

        List<Track> tracks = (new TrackLoader()).load(locator, genome);
        DataTrack igvTrack = (DataTrack) tracks.get(0);

        Chromosome chrom = handler.getChromosomeFromName("chr1");
        int res = 1;
        int binCount = (int) (chrom.getLength() / res + 1);
        int zoom = Math.max(0, (int) (Math.log(binCount / 700f) / log2));

        List<LocusScore> loci = getLocusScores(igvTrack, "chr1", 10000, 11000, zoom,
                WindowFunction.max);

        for (LocusScore locus : loci) {
            System.out.println(locus.getStart() + " " + locus.getEnd() + " " + locus.getScore());
        }
    }

    private static void aggregate(String inputBigWig, String bedFile, ChromosomeHandler handler, int resolution, String outFile, String genomeID) {
        List<Chromosome> chrList = new ArrayList<>(Arrays.asList(handler.getChromosomeArray()));
        List<org.broad.igv.feature.Chromosome> igvChrList = new ArrayList<>();
        for (Chromosome chrom : chrList) {
            igvChrList.add(chrom.toIGVChromosome());
        }

        Genome genome = new Genome(genomeID, igvChrList);
        ResourceLocator locator = new ResourceLocator(inputBigWig);

        List<Track> tracks = (new TrackLoader()).load(locator, genome);
        DataTrack igvTrack = (DataTrack) tracks.get(0);

        List<LocusScore> loci = getLocusScores(igvTrack, "chr1", 0, 100, 1000000,
                WindowFunction.mean);

        for (LocusScore locus : loci) {
            System.out.println(locus.getStart() + " " + locus.getEnd() + " " + locus.getScore());
        }


    }

    protected static List<LocusScore> getLocusScores(DataTrack igvTrack, String chr, long gStart, long gEnd, int zoom,
                                                     WindowFunction windowFunction) {
        igvTrack.setWindowFunction(windowFunction);
        org.broad.igv.track.LoadedDataInterval<List<LocusScore>> scores = igvTrack.getSummaryScores(chr, (int) gStart, (int) gEnd, zoom);
        // Problems with human not having the "chr".  Return scores if not 0, otherwise try adding "chr"
        if (scores.getFeatures().size() > 0) return scores.getFeatures();
        else return igvTrack.getSummaryScores("chr" + chr, (int) gStart, (int) gEnd, zoom).getFeatures();
    }

}
