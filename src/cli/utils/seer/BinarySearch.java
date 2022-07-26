package cli.utils.seer;

// Finds closest term in Binary search if the key we're searching for isn't in the sequence.
public class BinarySearch {
    public static int runBinarySearchIteratively(double[] sortedArray, double key, int low, int high) {
        //double diff = Integer.MAX_VALUE;
        if (key < sortedArray[low]) return low;
        if (key >= sortedArray[high]) return high;
        //int value = 0; use if you want value instead
        while (low <= high) {
            if (low == high) {
                return low;
            }
            int mid = low + ((high - low) / 2);
            // Checks distance of current term. Keeps track of the previous smallest distance from key.
            // If the distance is greater, updates distance and index.
            if (sortedArray[mid] < key && key <= sortedArray[mid + 1]) {
                //value = sortedArray[mid]; use this if you want the value
                return mid;
                // standard binary search operations to progress through the list.
            } else if (sortedArray[mid + 1] < key) {
                low = mid + 1;
            } else if (sortedArray[mid] >= key) {
                high = mid;
            }
        }
        // to ensure left_val < key <= right_value
        /* while (key > sortedArray[index]){
            index = index + 1;
        }
        */
        return -1;
    }

    public static void main(String[] args) {
        //Random rand = new Random();
        int[] numbers = new int[10000];
        int scalar = 20;
        // create random number generator w/ random values
        for (int i = 0; i < numbers.length; i++) {
            numbers[i] = (int) (scalar * Math.random());
        }
        double[] cumSums = SeerUtils.convertToCDF(numbers);

        System.out.println("Test 1: ");
        // test #1, testing w/ 100 random terms and over entire cdf
        // cases where test1 doesn't return true: target <= cumSums[index]
        for (int i = 0; i < 1000; i++) {
            // generate random num btwn 0-1, find closest index to that number
            double target = Math.random();
            int index = runBinarySearchIteratively(cumSums, target, 0, cumSums.length - 1);
            try {
                if (index >= 0 && index < cumSums.length) {
                    System.out.println("Target: " + target + " Index:" + cumSums[index] + " index + 1:" + cumSums[index + 1] +
                            "pass? " + (target <= cumSums[index + 1] && target > cumSums[index]));
                } else {
                    System.out.println("Target: " + target);
                }
            } catch (Exception e) {
                System.out.println("Target: " + target);
            }
        }

        System.out.println("Test 2: ");
        // test #2, testing w/ 100 random terms on a portion of cdf
        int randStartBin = (int) (Math.random() * cumSums.length);
        int randEndBin = randStartBin + (int) ((Math.random() * (cumSums.length - randStartBin)));
        double x0 = cumSums[randStartBin];
        double xF = cumSums[randEndBin];
        double range = xF - x0;
        for (int i = 0; i < 100; i++) {
            double r = Math.random();
            double target = r * range + x0;
            int index = runBinarySearchIteratively(cumSums, target, randStartBin, randEndBin);
            System.out.println("Target: " + target + "startBin: " + cumSums[randStartBin] + "endBin " +
                    cumSums[randEndBin] + " pass? " + (target <= cumSums[index + 1] && target > cumSums[index]));
        }
    }
}