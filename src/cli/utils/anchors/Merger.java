package cli.utils.anchors;

public class Merger {
    /*
    public static void callMergeAnchors(GenomeWideList<GenericLocus> locusList) {
        mergeAnchorsTakeSmaller(locusList);
    }


    private static void mergeAnchorsTakeSmaller(GenomeWideList<GenericLocus> anchorList) {
        anchorList.filterLists(new FeatureFilter<GenericLocus>() {
            @Override
            public List<GenericLocus> filter(String chr, List<GenericLocus> anchorList) {
                return mergeTakeSmaller(anchorList);
            }
        });
    }

    public static List<GenericLocus> mergeTakeSmaller(List<GenericLocus> anchors) {
        Collections.sort(anchors);

        Set<GenericLocus> merged = new HashSet<>();
        if (anchors.size() > 0) {
            GenericLocus current = (GenericLocus) anchors.get(0).deepClone();

            for (GenericLocus anchor : anchors) {
                if (anchor.hasOverlapWith(current)) {
                    current.mergeWithTakeSmaller(anchor);
                } else {
                    merged.add(current);
                    current = (GenericLocus) anchor.deepClone();
                }
            }
            merged.add(current); // in case last merger missed (i.e. boolean evaluated to true)
        }
        return new ArrayList<>(merged);
    }
     */
}
