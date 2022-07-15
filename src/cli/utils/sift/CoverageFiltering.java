package cli.utils.sift;

import javastraw.reader.block.ContactRecord;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class CoverageFiltering {
    public static void inPlaceFilterByNorms(Set<ContactRecord> initialPoints, double[] vector1, double[] vector1b, int scalar) {
        List<ContactRecord> toDelete = new LinkedList<>();
        for (ContactRecord record : initialPoints) {
            if (!normsAreGood(record, vector1, vector1b, scalar)) {
                toDelete.add(record);
            }
        }
        toDelete.forEach(initialPoints::remove);
    }

    public static boolean normsAreGood(ContactRecord record, double[] vector1, double[] vector2, int scalar) {
        return normIsGood(vector1, record, scalar) && normIsGood(vector2, record, scalar);
    }

    public static boolean normIsGood(double[] vector, ContactRecord record, int scalar) {
        return vector[record.getBinX() / scalar] > 1 && vector[record.getBinY() / scalar] > 1;
    }
}
