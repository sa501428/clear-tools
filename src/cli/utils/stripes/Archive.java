package cli.utils.stripes;

public class Archive {

    /*
    private boolean enrichedOEOverHorizontalWindow(int i, int j) {
        float[] rowOE = new float[3];
        Arrays.fill(rowOE, 1);
        for(int k = j - 1; k < j+2; k++){
            rowOE[0] *= getOEValue(i - 1, k);
            rowOE[1] *= getOEValue(i, k);
            rowOE[2] *= getOEValue(i + 1, k);
        }
        return rowOE[1] > 1.25 * rowOE[0] && rowOE[1] > 1.25 * rowOE[2];
    }

    private boolean enrichedOEOverVerticalWindow(int i, int j) {
        float[] colOE = new float[3];
        Arrays.fill(colOE, 1);
        for(int k = i - 1; k < i+2; k++){
            colOE[0] *= getOEValue(k, j-1);
            colOE[1] *= getOEValue(k, j);
            colOE[2] *= getOEValue(k, j+1);
        }
        return colOE[1] > 1.25 * colOE[0] && colOE[1] > 1.25 * colOE[2];
    }

    /*
    private void simpleHorizontalCall(int i, List<Feature2D> stripes) {
        int stripeStart = -1;
        int stripeEnd = -1;
        int stripeLen = 0;
        for (int j = i + minPeakDist; j < i + maxPeakDist; j++) {
            if (enrichedOverHorizontalWindow(i, j)) {
                if (stripeStart == -1) {
                    stripeStart = j;
                }
                stripeEnd = j;
                stripeLen++;
            } else {
                if (stripeLen >= minLengthStripe) {
                    stripes.add(makeHorizontalStripe(i, stripeStart, stripeEnd, resolution));
                }
                stripeStart = -1;
                stripeEnd = -1;
                stripeLen = 0;
            }
        }
        if (stripeLen >= minLengthStripe) {
            stripes.add(makeHorizontalStripe(i, stripeStart, stripeEnd, resolution));
        }
    }


    private void simpleVerticalCall(int j, List<Feature2D> stripes) {
        int stripeStart = -1;
        int stripeEnd = -1;
        int stripeLen = 0;
        for (int i = j - maxPeakDist; i < j - minPeakDist; i++) {
            if (enrichedOverVerticalWindow(i, j)) {
                if (stripeStart == -1) {
                    stripeStart = i;
                }
                stripeEnd = i;
                stripeLen++;
            } else {
                if (stripeLen >= minLengthStripe) {
                    stripes.add(makeVerticalStripe(j, stripeStart, stripeEnd, resolution));
                }
                stripeStart = -1;
                stripeEnd = -1;
                stripeLen = 0;
            }
        }
        if (stripeLen >= minLengthStripe) {
            stripes.add(makeVerticalStripe(j, stripeStart, stripeEnd, resolution));
        }
    }
    */
}
