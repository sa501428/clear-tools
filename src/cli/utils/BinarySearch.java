package cli.utils;

import java.util.Random;

// Finds closest term in Binary search if the key we're searching for isn't in the sequence.
public class BinarySearch {
    static public String runBinarySearchIteratively(double[] sortedArray, double key, int low, int high) {
        double diff = Integer.MAX_VALUE;
        //int value = 0; use if you want value instead
        int index = 0;
        while (low <= high) {
            int mid = low + ((high - low) / 2);
            // Checks distance of current term. Keeps track of the previous smallest distance from key.
            // If the distance is greater, updates distance and index.
            if (Math.abs(key - sortedArray[mid]) < diff) {
                diff = Math.abs(key - sortedArray[mid]);
                //value = sortedArray[mid]; use this if you want the value
                index = mid;
            }
            // standard binary search operations to progress through the list.
            if (sortedArray[mid] < key) {
                low = mid + 1;
            } else if (sortedArray[mid] > key){
                high = mid - 1;
            }
        }
        return "index - 1:" + (index - 1) + " index:" + (index) + " index + 1:" + (index + 1);
    }

    public static void main(String[] args) {
        Random rand = new Random();
        double[] randomNums = new double[100];
        double[] numSums = new double[100];

        // create random number generator w/ random values
        for (int i = 0; i < randomNums.length; i++) {
            randomNums[i] = rand.nextDouble();
        }
        // matrix of cumulative sums of the first array
        numSums[0] = randomNums[0];
        for (int i = 1; i < randomNums.length; i++) {
            numSums[i] = numSums[i - 1] + randomNums[i];
        }
        // normalizes every term by having last value = 1
        for (int i = 1; i < numSums.length; i++) {
            numSums[i] = numSums[i] / numSums[numSums.length - 1];
        }
        // generate random num btwn 0-1, find closest index to that number
        double target = Math.random();
        // print indexes close to the index we find (one below, one higher and index).
        System.out.println(runBinarySearchIteratively(numSums, target, 0, 100));
    }
}