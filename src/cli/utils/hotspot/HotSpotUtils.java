package cli.utils.hotspot;

import cli.utils.general.FusionTools;
import javastraw.reader.block.ContactRecord;

import java.util.*;

public class HotSpotUtils {
    public static void coalesceAndRetainCentroids(Set<ContactRecord> features,
                                                  int binRadius) {
        LinkedList<ContactRecord> featureLL = new LinkedList<>(features);
        List<ContactRecord> coalesced = new ArrayList<>();

        featureLL.sort((o1, o2) -> Float.compare(o1.getCounts(), o2.getCounts()));
        Collections.reverse(featureLL);

        while (!featureLL.isEmpty()) {
            ContactRecord pixel = featureLL.pollFirst();
            if (pixel != null) {
                coalesced.add(pixel);
                featureLL.remove(pixel);
                long x = pixel.getBinX();
                long y = pixel.getBinY();

                Set<ContactRecord> toRemove = new HashSet<>();
                for (ContactRecord px : featureLL) {
                    if (FusionTools.distance(x - px.getBinX(), y - px.getBinY()) <= binRadius) {
                        toRemove.add(px);
                    }
                }
                featureLL.removeAll(toRemove);
            }
        }

        features.clear();
        features.addAll(coalesced);
    }
}
