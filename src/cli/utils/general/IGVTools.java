package cli.utils.general;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.track.*;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IGVTools {
    private final static double log2 = Math.log(2.0D);

    public static List<LocusScore> getLocusScores(DataTrack igvTrack, Chromosome chromosome, long gStart, long gEnd, int resolution,
                                                  WindowFunction windowFunction) {
        int binCount = (int) (chromosome.getLength() / resolution + 1);
        int zoom = Math.max(0, (int) (Math.log(binCount / 700f) / log2));

        igvTrack.setWindowFunction(windowFunction);
        LoadedDataInterval<List<LocusScore>> scores = igvTrack.getSummaryScores(chromosome.getName(), (int) gStart, (int) gEnd, zoom);
        // Problems with human not having the "chr".  Return scores if not 0, otherwise try adding "chr"
        if (scores.getFeatures().size() > 0) return scores.getFeatures();
        else
            return igvTrack.getSummaryScores("chr" + chromosome.getName(), (int) gStart, (int) gEnd, zoom).getFeatures();
    }

    public static DataTrack loadBigWig(String inputBigWig, ChromosomeHandler handler) {
        Genome genome = getGenome(handler);
        ResourceLocator locator = new ResourceLocator(inputBigWig);
        List<Track> tracks = (new TrackLoader()).load(locator, genome);
        return (DataTrack) tracks.get(0);
    }

    public static TribbleFeatureSource loadBed(String inputBed, ChromosomeHandler handler) throws IOException, TribbleIndexNotFoundException {
        ResourceLocator locator = new ResourceLocator(inputBed);
        Genome genome = getGenome(handler);
        return TribbleFeatureSource.getFeatureSource(locator, genome);
    }

    static Genome getGenome(ChromosomeHandler handler) {
        List<Chromosome> chrList = new ArrayList<>(Arrays.asList(handler.getChromosomeArray()));
        List<org.broad.igv.feature.Chromosome> igvChrList = new ArrayList<>();
        for (Chromosome chrom : chrList) {
            igvChrList.add(toIGVChromosome(chrom));
        }
        return new Genome(handler.getGenomeID(), igvChrList);
    }

    public static org.broad.igv.feature.Chromosome toIGVChromosome(Chromosome chromosome) {
        return new org.broad.igv.feature.Chromosome(chromosome.getIndex(), chromosome.getName(),
                (int) chromosome.getLength()); // todo assumed for IGV
    }
}
