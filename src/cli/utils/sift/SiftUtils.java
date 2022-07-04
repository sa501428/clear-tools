package cli.utils.sift;

import cli.utils.FusionUtils;
import javastraw.reader.block.ContactRecord;

import java.util.*;

public class SiftUtils {
    public static List<ContactRecord> coalescePixelsToCentroid(Set<ContactRecord> features,
                                                               int resolution, int gRadius) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        LinkedList<ContactRecord> featureLL = new LinkedList<>(features);
        List<ContactRecord> coalesced = new ArrayList<>();

        featureLL.sort((o1, o2) -> Float.compare(o1.getCounts(), o2.getCounts()));
        Collections.reverse(featureLL);

        int clusterRadius = gRadius / resolution;

        while (!featureLL.isEmpty()) {
            ContactRecord pixel = featureLL.pollFirst();
            if (pixel != null) {
                coalesced.add(pixel);
                featureLL.remove(pixel);
                long x = pixel.getBinX();
                long y = pixel.getBinY();

                Set<ContactRecord> toRemove = new HashSet<>();
                for (ContactRecord px : featureLL) {
                    if (FusionUtils.distance(x - px.getBinX(),
                            y - px.getBinY()) <= clusterRadius) {
                        toRemove.add(px);
                    }
                }
                featureLL.removeAll(toRemove);
            }
        }

        return coalesced;
    }
}
