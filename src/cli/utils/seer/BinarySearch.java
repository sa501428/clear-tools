package cli.utils.seer;

import java.util.Random;

// Finds closest term in Binary search if the key we're searching for isn't in the sequence.
public class BinarySearch {
    static public String runBinarySearchIteratively(double[] sortedArray, double key, int low, int high) {
        double diff = Integer.MAX_VALUE;

        //if(key < sortedArray[0]) return 0;

        //int value = 0; use if you want value instead
        int index = 0;
        while (low <= high) {
            int mid = low + ((high - low) / 2);
            // Checks distance of current term. Keeps track of the previous smallest distance from key.
            // If the distance is greater, updates distance and index.

            // todo @Allen left_val < key <= right_value

            if (Math.abs(key - sortedArray[mid]) < diff) {
                diff = Math.abs(key - sortedArray[mid]);
                //value = sortedArray[mid]; use this if you want the value
                index = mid;
            }
            // standard binary search operations to progress through the list.
            if (sortedArray[mid] < key) {
                low = mid + 1;
            } else if (sortedArray[mid] > key) {
                high = mid - 1;
            }
        }
        return "Target: " + key + " index - 1:" + sortedArray[index - 1] + " index:" + sortedArray[index] +
                " pass? " + (key <= sortedArray[index] && key > sortedArray[index - 1]);
    }

    public static void main(String[] args) {
        Random rand = new Random();
        double[] numbers = new double[10000];
        double[] cumSums = new double[10000];

        int scalar = 20;

        // create random number generator w/ random values
        for (int i = 0; i < numbers.length; i++) {
            numbers[i] = scalar * rand.nextDouble();
        }
        // matrix of cumulative sums of the first array
        cumSums[0] = numbers[0];
        for (int i = 1; i < numbers.length; i++) {
            cumSums[i] = cumSums[i - 1] + numbers[i];
        }

        //
        double sum = cumSums[cumSums.length - 1];

        // normalizes every term by having last value = 1
        for (int i = 0; i < cumSums.length; i++) {
            cumSums[i] /= sum;
        }

        for (int i = 0; i < 100; i++) {
            // generate random num btwn 0-1, find closest index to that number
            double target = Math.random();
            // print indexes close to the index we find (one below, one higher and index).
            System.out.println(runBinarySearchIteratively(cumSums, target, 0, cumSums.length - 1));
        }
    }
}