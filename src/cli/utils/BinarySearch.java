package cli.utils;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Random;

// Finds closest term in Binary search if the key we're searching for isn't in the sequence.
public class BinarySearch {
    static public int runBinarySearchIteratively(int[] sortedArray, int key, int low, int high) {
        int diff = Integer.MAX_VALUE;
        //int value = 0; use if you want value instead
        int index = 0;
        while (low <= high) {
            int mid = low + ((high - low) / 2);
            // Checks distance of current term. Keeps track of the previous smallest distance from key.
            // If the distance is greater, updates distance and index.
            if (Math.abs(key-sortedArray[mid]) < diff) {
                diff = Math.abs(key-sortedArray[mid]);
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
        /*
        // in the case that the binary search doesn't find the term, index will still be max value
        if (index == Integer.MAX_VALUE){
            // find the value with minimum distance to key
            int newArray [] = new int[sortedArray.length];
            for(i = 0; i < sortedArray.length; i++){
                newArray[i] = Math.abs(sortedArray[i] - key);
            }
            // find the minimum difference. Not sure why this is giving an error as List is a type of collection?
            int min = Collections.min(Arrays.asList(newArray));
            low = Collections.min(Arrays.asList(sortedArray));
            high = Collections.max(Arrays.asList(sortedArray));

            // number we're looking for could either be distance above it or distance below (or both could be in).
            // I assume it doesn't matter which one we get.
            index = runBinarySearchIteratively(sortedArray, key + min low, high);
            // if index is still max value that means it wasn't distance above it. Check distance below
            if (index == Integer.MAX_VALUE){
                index = runBinarySearchIteratively(sortedArray, key - min, low, high);
            }
        }
        */
        return index;
    }

    public static void main(String[] args) {
        // create random number generator w/ random values btwn 0, 5 (a lot of values)
        // matrix of cumulative sums of the first array --> normalize by having last value = 1
        // generate random num btwn 0-1, find closest index to that number
        // print indexes close to the index we find (one below, one higher and index).
        System.out.println(runBinarySearchIteratively(new int[]{1,3,9,12,50,86,90}, 89, 0, 7));
    }
}