/**
 * 
 * Title: Program 4 
 * Files: ComparisonSort.java, Questions.txt;
 * 
 * Semester: Fall 2016 
 * 
 * Author: Katrina Van Laan 
 * Email:vanlaan@wisc.edu
 * Lecturer's Name: Charles Fischer
 *
 */

/**
 * This class implements six different comparison sorts (and an optional seventh
 * sort for extra credit):
 * <ul>
 * <li>selection sort</li>
 * <li>insertion sort</li>
 * <li>merge sort</li>
 * <li>quick sort</li>
 * <li>heap sort</li>
 * <li>selection2 sort</li>
 * <li>(extra credit) insertion2 sort</li>
 * </ul>
 * It also has a method that runs all the sorts on the same input array and
 * prints out statistics.
 */

public class ComparisonSort {

	// int which stores dataMoves performed on SortObjects
	private static int dataMoves = 0;

	/**
	 * Sorts the given array using the selection sort algorithm. You may use
	 * either the algorithm discussed in the on-line reading or the algorithm
	 * discussed in lecture (which does fewer data moves than the one from the
	 * on-line reading). Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 * @return
	 * @return
	 */

	/**
	 * Resets the datamoves counter to 0.
	 */
	private static void resetDataMoves() {
		dataMoves = 0;
	}

	/**
	 * Returns the number of datamoves that have been performed by a sort
	 * method.
	 *
	 * @return the number of datamoves that have performed
	 */
	public static int getDataMoves() {
		return dataMoves;
	}

	/**
	 * Swaps the elements in two given positions on a given array.
	 *
	 * @param A
	 *            the array the swap will be performed on
	 * @param i
	 *            the index of one element which will be swapped
	 * @param j
	 *            the index of another element which will be swapped
	 */
	private static <E extends Comparable<E>> void swap(E[] a, int i, int j) {
		E temp = a[i];
		a[i] = a[j];
		a[j] = temp;
		dataMoves = dataMoves + 3;
	}

	/**
	 * Sorts the given array using the selection sort algorithm. You may use
	 * either the algorithm discussed in the on-line reading or the algorithm
	 * discussed in lecture (which does fewer data moves than the one from the
	 * on-line reading). Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void selectionSort(E[] A) {

		for (int i = 0; i < A.length - 1; i++) {
			int minPos = i;
			for (int j = i + 1; j < A.length; j++) {
				if (A[j].compareTo(A[minPos]) < 1) {
					minPos = j;
				}
			}
			swap(A, minPos, i);
		}

	}

	/**
	 * Sorts the given array using the insertion sort algorithm. Note: after
	 * this method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void insertionSort(E[] A) {

		int k, j;
		E tmp;
		int N = A.length;

		for (k = 1; k < N; k++) {
			tmp = A[k];
			dataMoves++;
			j = k - 1;
			while ((j >= 0) && (A[j].compareTo(tmp) > 0)) {
				A[j + 1] = A[j]; // move one value over one place to the right
				dataMoves++;
				j--;
			}
			A[j + 1] = tmp; // insert kth value in correct place relative
			dataMoves++; // to previous values
		}
	}

	/**
	 * Sorts the given array using the merge sort algorithm. Note: after this
	 * method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void mergeSort(E[] A) {

		mergeAux(A, 0, A.length - 1);
	}

	/**
	 * Helper method for mergesort
	 * 
	 * @param A
	 *            array to be sorted
	 * @param low
	 *            low index
	 * 
	 * @param high
	 *            high index
	 */
	private static <E extends Comparable<E>> void mergeAux(E[] A, int low, int high) {

		if (low == high)
			return;

		int mid = (low + high) / 2;
		mergeAux(A, low, mid);
		mergeAux(A, mid + 1, high);
		E[] tmp = (E[]) (new Comparable[high - low + 1]);
		dataMoves++;
		int left = low;
		int right = mid + 1;
		int pos = 0;
		while ((left <= mid) && (right <= high)) {

			if (A[left].compareTo(A[right]) <= 0) {
				tmp[pos] = A[left];
				dataMoves++;
				left++;
			} else {
				tmp[pos] = A[right];
				dataMoves++;
				right++;
			}
			pos++;
		}
		while (left <= mid) {
			tmp[pos] = A[left];
			dataMoves++;
			left++;
			pos++;
		}
		while (right <= high) {
			tmp[pos] = A[right];
			dataMoves++;
			right++;
			pos++;
		}
		System.arraycopy(tmp, 0, A, low, tmp.length);
		dataMoves = dataMoves + tmp.length;
	}

	/**
	 * Sorts the given array using the quick sort algorithm, using the median of
	 * the first, last, and middle values in each segment of the array as the
	 * pivot value. Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */

	public static <E extends Comparable<E>> void quickSort(E[] A) {
		if (A == null || A.length == 0) {
			return;
		}
		int length = A.length;
		quickSortHelp(A, 0, length - 1);
	}

	/**
	 * Helper Function to recursively implement quicksort
	 * 
	 * @param input
	 *            the array to be sorted
	 * @param low
	 *            to be used as lower index
	 * 
	 * @param high
	 *            to be used as higher index
	 */
	private static <E extends Comparable<E>> void quickSortHelp(E[] input, int low, int high) {
		int i = low;
		int j = high;
		E pivot = input[low + (high - low) / 2];
		dataMoves++;
		// Divide into two arrays
		while (i <= j) {
			while (input[i].compareTo(pivot) < 0) {
				i++;
			}
			while (input[j].compareTo(pivot) > 0) {

				j--;
			}
			if (i <= j) {
				swap(input, i, j); // move index to next position on both sides
				i++;
				//
				j--;
			}
		} // calls quickSort() method recursively
		if (low < j) {
			quickSortHelp(input, low, j);
		}
		if (i < high) {
			quickSortHelp(input, i, high);
		}
	}

	/**
	 * Sorts the given array using the heap sort algorithm outlined below. Note:
	 * after this method finishes the array is in sorted order.
	 * <p>
	 * The heap sort algorithm is:
	 * </p>
	 * 
	 * <pre>
	 * for each i from 1 to the end of the array
	 *     insert A[i] into the heap (contained in A[0]...A[i-1])
	 *     
	 * for each i from the end of the array up to 1
	 *     remove the max element from the heap and put it in A[i]
	 * </pre>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */

	public static <E extends Comparable<E>> void heapSort(E[] A) {
		heapify(A);
		for (int i = N; i > 0; i--) {
			swap(A, 0, i);
			N = N - 1;
			maxheap(A, 0);
		}
	}

	// helper field for heapSort
	private static int N;

	/**
	 * Helper Function to build a heap
	 * 
	 * @param arr
	 *            the array the heap will be built off of;
	 */
	public static <E extends Comparable<E>> void heapify(E[] arr) {
		N = arr.length - 1;
		for (int i = N / 2; i >= 0; i--) {
			maxheap(arr, i);
		}
	}

	/**
	 * Helper Function to swap largest element in heap
	 * 
	 * @param E
	 *            arr[] the array where swapping will occur;
	 * @param i
	 *            index which will be swapped with largest
	 */
	public static <E extends Comparable<E>> void maxheap(E arr[], int i) {
		int left = 2 * i;
		int right = 2 * i + 1;
		int max = i;
		if (left <= N && arr[left].compareTo(arr[i]) > 0)
			max = left;
		if (right <= N && arr[right].compareTo(arr[max]) > 0)
			max = right;

		if (max != i) {
			swap(arr, i, max);
			maxheap(arr, max);
		}
	}

	/**
	 * Sorts the given array using the selection2 sort algorithm outlined below.
	 * Note: after this method finishes the array is in sorted order.
	 * <p>
	 * The selection2 sort is a bi-directional selection sort that sorts the
	 * array from the two ends towards the center. The selection2 sort algorithm
	 * is:
	 * </p>
	 * 
	 * <pre>
	 * begin = 0, end = A.length-1
	 * 
	 * // At the beginning of every iteration of this loop, we know that the 
	 * // elements in A are in their final sorted positions from A[0] to A[begin-1]
	 * // and from A[end+1] to the end of A.  That means that A[begin] to A[end] are
	 * // still to be sorted.
	 * do
	 *     use the MinMax algorithm (described below) to find the minimum and maximum 
	 *     values between A[begin] and A[end]
	 *     
	 *     swap the maximum value and A[end]
	 *     swap the minimum value and A[begin]
	 *     
	 *     ++begin, --end
	 * until the middle of the array is reached
	 * </pre>
	 * <p>
	 * The MinMax algorithm allows you to find the minimum and maximum of N
	 * elements in 3N/2 comparisons (instead of 2N comparisons). The way to do
	 * this is to keep the current min and max; then
	 * </p>
	 * <ul>
	 * <li>take two more elements and compare them against each other</li>
	 * <li>compare the current max and the larger of the two elements</li>
	 * <li>compare the current min and the smaller of the two elements</li>
	 * </ul>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */

	public static <E extends Comparable<E>> void selection2Sort(E[] A) {
		int begin = -1;
		int end = A.length - 2;
		boolean swapped;
		do {
			swapped = false;
			for (int i = ++begin; i <= end; i++) {
				if (A[i].compareTo(A[i + 1]) > 0) {
					swap(A, i, i + 1);
					swapped = true;
				}
			}

			if (!swapped) {
				break;
			}

			swapped = false;
			for (int i = --end; i >= begin; i--) {
				if (A[i].compareTo(A[i + 1]) > 0) {
					swap(A, i, i + 1);
					swapped = true;
				}
			}
		} while (swapped);
	}

	/**
	 * <b>Extra Credit:</b> Sorts the given array using the insertion2 sort
	 * algorithm outlined below. Note: after this method finishes the array is
	 * in sorted order.
	 * <p>
	 * The insertion2 sort is a bi-directional insertion sort that sorts the
	 * array from the center out towards the ends. The insertion2 sort algorithm
	 * is:
	 * </p>
	 * 
	 * <pre>
	 * precondition: A has an even length
	 * left = element immediately to the left of the center of A
	 * right = element immediately to the right of the center of A
	 * if A[left] > A[right]
	 *     swap A[left] and A[right]
	 * left--, right++ 
	 *  
	 * // At the beginning of every iteration of this loop, we know that the elements
	 * // in A from A[left+1] to A[right-1] are in relative sorted order.
	 * do
	 *     if (A[left] > A[right])
	 *         swap A[left] and A[right]
	 *  
	 *     starting with with A[right] and moving to the left, use insertion sort 
	 *     algorithm to insert the element at A[right] into the correct location 
	 *     between A[left+1] and A[right-1]
	 *     
	 *     starting with A[left] and moving to the right, use the insertion sort 
	 *     algorithm to insert the element at A[left] into the correct location 
	 *     between A[left+1] and A[right-1]
	 *  
	 *     left--, right++
	 * until left has gone off the left edge of A and right has gone off the right 
	 *       edge of A
	 * </pre>
	 * <p>
	 * This sorting algorithm described above only works on arrays of even
	 * length. If the array passed in as a parameter is not even, the method
	 * throws an IllegalArgumentException
	 * </p>
	 *
	 * @param A
	 *            the array to sort
	 * @throws IllegalArgumentException
	 *             if the length or A is not even
	 */
	public static <E extends Comparable<E>> void insertion2Sort(E[] A) throws IllegalArgumentException {
		if (A.length % 2 != 0) {
			throw new IllegalArgumentException();
		}

		E leftTmp;
		E rightTmp;
		int length = A.length;
		int j = 0;
		int left = (length / 2) - 1;
		int right = (length / 2) + 1;
		int N = length / 2;

		for (int i = 1; i < N; i++) {
			for (int k = right; k > left; k--) {
				if (A[left].compareTo(A[right]) > 0) {
					swap(A, left, right);
				}
				rightTmp = A[k];
				dataMoves++;
				j = k - 1;
				while ((j > left) && (A[j].compareTo(rightTmp) > 0)) {
					A[j + 1] = A[j];
					j--;
					dataMoves++;
				}
				A[j + 1] = rightTmp;
				dataMoves++;
			}

			for (int k = left; k < right; k++) {
				leftTmp = A[k];
				dataMoves++;
				j = k + 1;
				while ((j < right) && (A[j].compareTo(leftTmp) < 0)) {
					A[j - 1] = A[j];
					j++;
					dataMoves++;
				}
				A[j - 1] = leftTmp;
				dataMoves++;
			}
			left--;
			right++;
		}
	}

	/**
	 * Internal helper for printing rows of the output table.
	 * 
	 * @param sort
	 *            name of the sorting algorithm
	 * @param compares
	 *            number of comparisons performed during sort
	 * @param moves
	 *            number of data moves performed during sort
	 * @param milliseconds
	 *            time taken to sort, in milliseconds
	 */
	private static void printStatistics(String sort, int compares, int moves, long milliseconds) {
		System.out.format("%-23s%,15d%,15d%,15d\n", sort, compares, moves, milliseconds);
	}

	/**
	 * Sorts the given array using the six (seven with the extra credit)
	 * different sorting algorithms and prints out statistics. The sorts
	 * performed are:
	 * <ul>
	 * <li>selection sort</li>
	 * <li>insertion sort</li>
	 * <li>merge sort</li>
	 * <li>quick sort</li>
	 * <li>heap sort</li>
	 * <li>selection2 sort</li>
	 * <li>(extra credit) insertion2 sort</li>
	 * </ul>
	 * <p>
	 * The statistics displayed for each sort are: number of comparisons, number
	 * of data moves, and time (in milliseconds).
	 * </p>
	 * <p>
	 * Note: each sort is given the same array (i.e., in the original order) and
	 * the input array A is not changed by this method.
	 * </p>
	 * 
	 * @param A
	 *            the array to sort
	 */
	static public void runAllSorts(SortObject[] A) {
		System.out.format("%-23s%15s%15s%15s\n", "algorithm", "data compares", "data moves", "milliseconds");
		System.out.format("%-23s%15s%15s%15s\n", "---------", "-------------", "----------", "------------");

		// TODO: run each sort and print statistics about what it did

		SortObject[] test1 = new SortObject[A.length];
		System.arraycopy(A, 0, test1, 0, A.length);

		int begTime = (int) System.currentTimeMillis();
		selectionSort(test1);
		int endTime = (int) System.currentTimeMillis();

		System.out.format("%-23s%15s%15s%15s\n", "selection", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);

		SortObject.resetCompares();
		ComparisonSort.resetDataMoves();

		SortObject[] test2 = new SortObject[A.length];
		System.arraycopy(A, 0, test2, 0, A.length);

		begTime = (int) System.currentTimeMillis();
		insertionSort(test2);
		endTime = (int) System.currentTimeMillis();
		System.out.format("%-23s%15s%15s%15s\n", "insertion", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);
		SortObject.resetCompares();

		ComparisonSort.resetDataMoves();

		SortObject[] test3 = new SortObject[A.length];
		System.arraycopy(A, 0, test3, 0, A.length);

		begTime = (int) System.currentTimeMillis();
		mergeSort(test3);
		endTime = (int) System.currentTimeMillis();
		System.out.format("%-23s%15s%15s%15s\n", "merge", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);
		SortObject.resetCompares();

		ComparisonSort.resetDataMoves();

		SortObject[] test4 = new SortObject[A.length];
		System.arraycopy(A, 0, test4, 0, A.length);

		begTime = (int) System.currentTimeMillis();
		quickSort(test4);
		endTime = (int) System.currentTimeMillis();

		System.out.format("%-23s%15s%15s%15s\n", "quick", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);

		ComparisonSort.resetDataMoves();
		SortObject.resetCompares();

		SortObject[] test5 = new SortObject[A.length];
		System.arraycopy(A, 0, test5, 0, A.length);

		begTime = (int) System.currentTimeMillis();
		heapSort(test5);
		endTime = (int) System.currentTimeMillis();

		System.out.format("%-23s%15s%15s%15s\n", "heap", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);

		ComparisonSort.resetDataMoves();

		SortObject[] test6 = new SortObject[A.length];
		System.arraycopy(A, 0, test6, 0, A.length);

		SortObject.resetCompares();
		begTime = (int) System.currentTimeMillis();
		selection2Sort(test6);
		endTime = (int) System.currentTimeMillis();

		System.out.format("%-23s%15s%15s%15s\n", "selection2Sort", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);

		SortObject.resetCompares();

		ComparisonSort.resetDataMoves();

		SortObject[] test7 = new SortObject[A.length];
		System.arraycopy(A, 0, test7, 0, A.length);

		begTime = (int) System.currentTimeMillis();
		insertion2Sort(test7);
		endTime = (int) System.currentTimeMillis();
		System.out.format("%-23s%15s%15s%15s\n", "insertion2Sort", String.format("%,d", SortObject.getCompares()),
				String.format("%,d", ComparisonSort.getDataMoves()), endTime - begTime);

		SortObject.resetCompares();

		ComparisonSort.resetDataMoves();
	}
}