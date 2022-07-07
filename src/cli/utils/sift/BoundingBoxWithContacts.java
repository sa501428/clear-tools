package cli.utils.sift;

import cli.utils.flags.Utils;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class BoundingBoxWithContacts {

    private final static int buffer = 10;
    private final static int width = 2;
    private final List<ContactRecord> contacts;
    private final int scalar;
    private int minR, maxR, minC, maxC;

    public BoundingBoxWithContacts(List<ContactRecord> contacts, int scalar) {
        this.contacts = contacts;
        this.scalar = scalar;
        minR = contacts.get(0).getBinX() / scalar;
        minC = contacts.get(0).getBinY() / scalar;
        maxR = minR;
        maxC = minC;
        setBounds();
    }

    private void setBounds() {
        for (ContactRecord contact : contacts) {
            minR = Math.min(minR, contact.getBinX() / scalar - buffer);
            minC = Math.min(minC, contact.getBinY() / scalar - buffer);

            maxR = Math.max(maxR, contact.getBinX() / scalar + buffer);
            maxC = Math.max(maxC, contact.getBinY() / scalar + buffer);
        }
        if (minR < 0) minR = 0;
        if (minC < 0) minC = 0;
    }

    public float[][] getRegionSpanned(MatrixZoomData zd, NormalizationType norm) {
        return Utils.getRegion(zd, minR, minC, maxR, maxC, norm);
    }

    public Set<ContactRecord> findPointsNotEnriched(MatrixZoomData zdLow, NormalizationType norm) {
        Set<ContactRecord> toRemove = new HashSet<>();
        float[][] region = getRegionSpanned(zdLow, norm);
        for (ContactRecord contact : contacts) {
            if (!isProperlyEnriched(region, contact.getBinX() / scalar - minR, contact.getBinY() / scalar - minC)) {
                toRemove.add(contact);
            }
        }
        return toRemove;
    }

    private boolean isProperlyEnriched(float[][] region, int i, int j) {
        float loopVal = region[i][j];
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int r = i - width; r <= i + width; r++) {
            for (int c = j - width; c <= j + width; c++) {
                if (r == i && c == j) continue;
                // too much "sparsity" near this region
                if (region[r][c] < 1e-10) return false;
                // something else nearby is bigger
                if (region[r][c] > loopVal) return false;

                stats.addValue(region[r][c]);
            }
        }

        // neighborhood stats
        double mean = stats.getMean();
        //double max = stats.getMax();
        double median = stats.getPercentile(50);

        return reasonablyEnrichedRelativeTo(mean, loopVal) &&
                //reasonablyEnrichedRelativeTo(max, loopVal) &&
                reasonablyEnrichedRelativeTo(median, loopVal);
    }

    private boolean reasonablyEnrichedRelativeTo(double denom, float num) {
        double val = num / denom;
        return val > 1.2 && val < 50;
    }
}
