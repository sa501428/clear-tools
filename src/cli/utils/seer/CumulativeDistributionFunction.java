package cli.utils.seer;

import cli.utils.sift.SimpleLocation;
import javastraw.expected.ExpectedUtils;
import javastraw.reader.block.ContactRecord;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

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
        toSave.clear();
    }

    private List<ContactRecord> extractTheRecordsWeWantToSave(Iterator<ContactRecord> normalizedIterator, int maxDist) {
        List<ContactRecord> toSave = new LinkedList<>();
        // get only counts within a certain range
        while (normalizedIterator.hasNext()) {
            ContactRecord record = normalizedIterator.next();
            if (record.getCounts() > 0 & ExpectedUtils.getDist(record) < maxDist) {
                toSave.add(record);
            }
        }
        return toSave;
    }

    private void populateCDFandLocations(double[] cdf, SimpleLocation[] genomeLocations, List<ContactRecord> toSave) {
        double total = 0;
        int index = 0;

        // fills in the cdf and genome locations
        for (ContactRecord record : toSave) {
            total += record.getCounts();
            cdf[index] = total;
            genomeLocations[index] = new SimpleLocation(record, resolution);
            index++;
        }

        // normalizes cdf
        for (int i = 0; i < cdf.length; i++) {
            cdf[i] /= total;
        }
    }

    public SimpleLocation createRandomPoint(Random rand) {
        double target = rand.nextDouble();
        int index = BinarySearch.runBinarySearchIteratively(cdf, target, 0, cdf.length - 1);
        return genomeLocations[index];
    }
}
