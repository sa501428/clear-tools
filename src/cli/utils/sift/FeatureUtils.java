package cli.utils.sift;

import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class FeatureUtils {

    public static List<Feature2D> convertToFeature2Ds(Set<ContactRecordBox> records, Chromosome c1) {
        List<Feature2D> features = new ArrayList<>();
        for (ContactRecordBox record : records) {
            features.add(record.toFeature2D(c1));
        }
        return features;
    }
}
