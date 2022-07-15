package cli.utils.seer;

import cli.utils.sift.SimpleLocation;
import javastraw.reader.block.ContactRecord;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class CumulativeDistributionFunction {

    private final int resolution;
    private final double[] cdf;
    private final SimpleLocation[] genomeLocations;

    public CumulativeDistributionFunction(Iterator<ContactRecord> normalizedIterator, int maxGenomeDist, int resolution) {
        this.resolution = resolution;
        List<ContactRecord> toSave = extractTheRecordsWeWantToSave(normalizedIterator, maxGenomeDist / resolution);
        cdf = new double[toSave.size()];
        genomeLocations = new SimpleLocation[toSave.size()];
        populateCDFandLocations(cdf, genomeLocations, toSave);
    }

    private List<ContactRecord> extractTheRecordsWeWantToSave(Iterator<ContactRecord> normalizedIterator, int maxDist) {
        List<ContactRecord> toSave = new LinkedList<>();

        // todo @Allen iterate on all the contacts
        // check that the counts > 0, and that the dist < maxDist
        // save ContactRecord to a valid list
        return toSave;
    }

    private void populateCDFandLocations(double[] cdf, SimpleLocation[] genomeLocations, List<ContactRecord> toSave) {
        double total = 0;
        int index = 0;

        for (ContactRecord record : toSave) {
            total += record.getCounts();
            cdf[index] = total;
            genomeLocations[index] = new SimpleLocation(record, resolution);
            index++;
        }

        for (int i = 0; i < cdf.length; i++) {
            cdf[i] /= total;
        }
    }
}
