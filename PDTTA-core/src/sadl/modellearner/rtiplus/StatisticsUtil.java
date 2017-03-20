/**
 * This file is part of SADL, a library for learning all sorts of (timed) automata and performing sequence-based anomaly detection.
 * Copyright (C) 2013-2016  the original author or authors.
 *
 * SADL is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * SADL is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with SADL.  If not, see <http://www.gnu.org/licenses/>.
 */
package sadl.modellearner.rtiplus;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.IntStream;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.Precision;

import gnu.trove.iterator.TDoubleIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

/**
 * A utility class that provides some common statistical functions.
 * 
 * @author Fabian Witter
 */
public abstract class StatisticsUtil {

	/**
	 * Returns the mean of a collection of {@link Number} objects.
	 * 
	 * @param values
	 *            {@link Collection} of values; elements may not be {@code null}; {@code NaN} and infinite values will be ignored
	 * 
	 * @return The mean of the given {@link Collection}
	 */
	public static double calculateMean(TDoubleList values) {

		int count = 0;
		double total = 0.0;
		final TDoubleIterator iterator = values.iterator();
		while (iterator.hasNext()) {
			final double value = iterator.next();
			if (!Double.isNaN(value) && !Double.isInfinite(value)) {
				total += value;
				count++;
			}
		}
		return total / count;
	}

	/**
	 * 0: Minimum regular value<br>
	 * 1: Q1<br>
	 * 2: Median<br>
	 * 3: Q3<br>
	 * 4: Maximum regular value<br>
	 * 5: Mean<br>
	 * 
	 * @param values
	 * @param copyAndSort
	 * @return
	 */
	public static double[] calculateBox(TDoubleList values, boolean copyAndSort) {

		final double[] box = new double[6];
		if (values != null && copyAndSort) {
			final TDoubleList copy = new TDoubleArrayList(values);
			copy.sort();
			values = copy;
		} else {
			throw new IllegalStateException("null not allowed as list for calculateBox");
		}
		box[5] = calculateMean(values);
		box[2] = calculateMedian(values, false);
		box[1] = calculateQ1(values, false);
		box[3] = calculateQ3(values, false);
		final double iqr = box[3] - box[1];
		box[0] = box[1] - (1.5 * iqr);
		box[4] = box[3] + (1.5 * iqr);
		return box;
	}

	public static double calculateQ1(TDoubleList values, boolean copyAndSort) {

		double result = Double.NaN;
		if (values != null) {
			if (copyAndSort) {
				final TDoubleList copy = new TDoubleArrayList(values);
				copy.sort();
				values = copy;
			}
			final int count = values.size();
			if (count > 0) {
				if (count % 2 == 1) {
					if (count > 1) {
						result = calcMedian(values, 0, (count - 3) / 2);
					} else {
						result = calcMedian(values, 0, 0);
					}
				} else {
					result = calcMedian(values, 0, (count / 2) - 1);
				}

			}
		}
		return result;
	}

	public static double calculateQ3(TDoubleList values, boolean copyAndSort) {

		double result = Double.NaN;
		if (values != null) {
			if (copyAndSort) {
				final TDoubleList copy = new TDoubleArrayList(values);
				copy.sort();
				values = copy;
			}
			final int count = values.size();
			if (count > 0) {
				if (count % 2 == 1) {
					if (count > 1) {
						result = calcMedian(values, (count + 1) / 2, count - 1);
					} else {
						result = calcMedian(values, 0, 0);
					}
				} else {
					result = calcMedian(values, count / 2, count - 1);
				}
			}
		}
		return result;
	}

	public static double calculateMedcoupleSlow(TDoubleList values, double median, boolean copyAndSort) {

		if (copyAndSort) {
			final TDoubleList copy = new TDoubleArrayList(values);
			copy.sort();
			copy.reverse();
			values = copy;
		}

		// values[0] contains max value
		final double scale = 2 * Math.abs(values.get(0));

		final TDoubleList vPlus = new TDoubleArrayList(values.size());
		final TDoubleList vMinus = new TDoubleArrayList(values.size());
		values.forEach(x -> {
			if (Precision.equals(x, median)) {
				vMinus.add((x - median) / scale);
				vPlus.add((x - median) / scale);
			} else if (x < median) {
				vMinus.add((x - median) / scale);
			} else {
				vPlus.add((x - median) / scale);
			}
			return true;
		});

		final int p = vPlus.size();
		final int q = vMinus.size();

		final BiFunction<Integer, Integer, Double> h = (i, j) -> {
			final double a = vPlus.get(i.intValue());
			final double b = vMinus.get(j.intValue());

			if (Precision.equals(a, b)) {
				return new Double(Math.signum(p - 1 - i.intValue() - j.intValue()));
			} else {
				return new Double((a + b) / (a - b));
			}
		};

		final TDoubleList H = new TDoubleArrayList(p * q);
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < q; j++) {
				H.add(h.apply(i, j));
			}
		}

		final double medcouple = calculateMedian(H, false);
		return medcouple;
	}

	public static double calculateMedcouple(TDoubleList values, double median, boolean copyAndSort) {

		System.out.println("start: " + values);

		if (copyAndSort) {
			final TDoubleList copy = new TDoubleArrayList(values);
			copy.sort();
			copy.reverse();
			values = copy;
		}

		// values[0] contains max value
		final double scale = 2 * Math.abs(values.get(0));

		final TDoubleList vPlus = new TDoubleArrayList(values.size());
		final TDoubleList vMinus = new TDoubleArrayList(values.size());
		values.forEach(x -> {
			if (Precision.equals(x, median)) {
				vMinus.add((x - median) / scale);
				vPlus.add((x - median) / scale);
			} else if (x < median) {
				vMinus.add((x - median) / scale);
			} else {
				vPlus.add((x - median) / scale);
			}
			return true;
		});

		final int p = vPlus.size();
		final int q = vMinus.size();

		final BiFunction<Integer, Integer, Double> h = (i, j) -> {
			final double a = vPlus.get(i.intValue());
			final double b = vMinus.get(j.intValue());

			if (Precision.equals(a, b)) {
				return new Double(Math.signum(p - 1 - i.intValue() - j.intValue()));
			} else {
				return new Double((a + b) / (a - b));
			}
		};

		System.out.println("###### H ####");
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < q; j++) {
				System.out.print(h.apply(i, j) + "\t");
			}
			System.out.println();
		}
		System.out.println("#############");

		// begin Kth pair algorithm (Johnson & Mizoguchi)

		// the initial left and right boundaries, two vectors of size p
		int[] l = new int[p];
		Arrays.fill(l, 0);
		int[] r = new int[p];
		Arrays.fill(r, q - 1);

		// number of entries to the left of the left boundary
		int lTotal = 0;
		// number of entries to the left of the right boundary
		int rTotal = p * q;

		// since we are indexing from zero, the medcouple index is one
		// less than its rank
		final int medcoupleIndex = (int) Math.floor(rTotal / 2.0);

		// iterate while the number of entries between the boundaries is
		// greater than the number of rows in the matrix
		while ((rTotal - lTotal) > p) {

			// compute row medians and their associated weights, but skip
			// any rows that are already empty
			final TIntList middleIdx = new TIntArrayList(p);
			for (int i = 0; i < l.length; i++) {
				if (l[i] <= r[i]) {
					// Row is not empty -> Add index
					middleIdx.add(i);
				}
			}
			final int[] finalL = l;
			final int[] finalR = r;
			final double[] rowMedians = Arrays.stream(middleIdx.toArray()).mapToDouble(i -> h.apply(i, (int) Math.floor((finalL[i] + finalR[i]) / 2.0)))
					.toArray();
			final double[] weights = Arrays.stream(middleIdx.toArray()).mapToDouble(i -> finalR[i] - finalL[i] + 1).toArray();

			final double WM = weightedMedian(rowMedians, weights);

			// new tentative right and left boundaries
			final int[] P = greaterThanH(h, p, q, WM);
			final int[] Q = lessThanH(h, p, q, WM);

			final int pTotal = IntStream.of(P).sum() + P.length;
			final int qTotal = IntStream.of(Q).sum();

			// determine which entries to discard, or if we've found the medcouple
			if (medcoupleIndex <= pTotal - 1) {
				r = P;
				rTotal = pTotal;
			} else if (medcoupleIndex > qTotal - 1) {
				l = Q;
				lTotal = qTotal;
			} else {
				// found the medcouple, rank of the weighted median equals medcouple index
				return WM;
			}
		}

		// did not find the medcouple, but there are very few tentative entries remaining
		final TDoubleList remaining = new TDoubleArrayList(p * q);
		for (int i = 0; i < p; i++) {
			for (int j = l[i]; j <= r[i]; j++) {
				remaining.add(h.apply(i, j));
			}
		}

		System.out.println(remaining);
		// select the medcouple by rank amongst the remaining entries
		final double medcouple = remaining.get(medcoupleIndex - lTotal);
		return medcouple;
	}

	public static double weightedMedian(double[] values, double[] weights) {
		return wquickSelect(Arrays.copyOf(values, values.length), Arrays.copyOf(weights, weights.length), values.length, -1.0, 0.0, 0.0);
	}

	/*
	 * This Quickselect routine is based on the algorithm described in
	 * "Numerical recipes in C", Second Edition,
	 * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	 * This code by Nicolas Devillard - 1998. Public domain.
	 */
	private static double wquickSelect(double[] x, double[] w, int N, double W0, double Wl, double Wr) {

		int low = 0;
		int high = N - 1;
		int middle, ll, hh;
		double wSumLeft = 0.0;
		double wLeft = Wl; // The sum of weights of discarded element to the left
		double wRight = Wr; // ... to the right

		if (W0 < 0.0) {
			// W0 is not provided -> compute
			W0 = 0.0;
			for (int i = 0; i < N; i++) {
				W0 += w[i];
			}
			W0 *= 0.5;
		}

		for (;;) {
			if (high <= low) { /* One element only */
				return x[low];
			}

			if (high == low + 1) { /* Two elements only */
				if (x[low] > x[high]) {
					swap(x, low, high);
					swap(w, low, high);
				}
				if (wLeft + w[low] >= W0) {
					return x[low];
				} else {
					return x[high];
				}
			}

			/* Find median of low, middle and high items; swap into position low */
			middle = (low + high) / 2;
			// Take 3-sample median of low,middle & high and put it in x[low]
			if (x[middle] > x[high]) {
				swap(x, middle, high);
				swap(w, middle, high);
			}
			if (x[low] > x[high]) {
				swap(x, low, high);
				swap(w, low, high);
			}
			// Now the largest sample is in x[high]
			if (x[middle] > x[low]) {
				swap(x, middle, low);
				swap(w, middle, low);
			}
			// x[middle] < x[low] < x[high]

			/* Swap low item (now in position middle) into position (low+1) */
			swap(x, middle, low + 1);
			swap(w, middle, low + 1);

			// We will use x[low] as the pivot
			final double piv = x[low];

			/* Nibble from each end towards middle, swapping items when stuck */
			ll = low + 1;
			hh = high;
			for (;;) {
				do {
					wSumLeft += w[ll];
					ll++;
				} while (piv > x[ll]);
				do {
					hh--;
				} while (x[hh] > piv);

				if (hh < ll) {
					break;
				}

				// Swap those items and their corresponding weight
				swap(x, ll, hh);
				swap(w, ll, hh);
			}

			/* Swap middle item (in position low) back into correct position */
			swap(x, low, hh);
			swap(w, low, hh);
			// wSumLeft += w[hh]; Do not include the pivot

			// In case there were elements equal to the pivot, we need to adjust the sum of weights
			while (ll > hh + 1) {
				ll--;
				wSumLeft -= w[ll];
			}

			// Note: wSumLeft is the sum of all weights left of the pivot excluding the
			// pivot itself.
			/* Re-set active partition */
			if ((wLeft + wSumLeft) <= W0) {
				if ((wLeft + wSumLeft + w[hh]) > W0) {
					// The pivot is the WM -> return
					return x[hh];
				}
				// The pivot was smaller than the WM -> discard elements to the left
				low = hh + 1; // Set the left limit right to the pivot
				wLeft += w[hh] + wSumLeft;
			} else {
				// The pivot was larger than the WM -> discard elements to the right
				high = hh - 1; // Set the right limit left to the pivot
				wRight = 2 * W0 - wSumLeft - wLeft;
			}
			// Reset the counter
			wSumLeft = 0.0;
		}
	}

	private static void swap(double[] x, int a, int b) {
		final double temp = x[a];
		x[a] = x[b];
		x[b] = temp;
	}

	public static double weightedMedianSlow(double[] values, double[] weights) {

		if (values.length != weights.length) {
			throw new IllegalArgumentException("Values and weigths must have the same size");
		}

		if (Arrays.stream(weights).anyMatch(w -> w < 0.0)) {
			throw new IllegalArgumentException("All weights must be at least 0");
		}

		final double weightSum = Arrays.stream(weights).sum();

		final List<Pair<Double, Double>> weightedVals = new ArrayList<>(values.length);
		for (int i = 0; i < values.length; i++) {
			weightedVals.add(Pair.of(new Double(values[i]), new Double(weights[i])));
		}

		Collections.sort(weightedVals, (p1, p2) -> p1.getLeft().compareTo(p2.getLeft()));

		int k = 0;
		double sum = weightSum - weightedVals.get(0).getRight().doubleValue();
		// sum is the total weight of all `x[i] > x[k]`

		while (sum > weightSum / 2.0) {
			++k;
			sum -= weightedVals.get(k).getRight().doubleValue();
		}

		return weightedVals.get(k).getLeft().doubleValue();
	}

	private static int[] greaterThanH(BiFunction<Integer, Integer, Double> h, int p, int q, double u) {
		// h is the kernel function, h(i,j) gives the ith, jth entry of H
		// p and q are the number of rows and columns of the kernel matrix H

		// vector of size p
		final int[] P = new int[p];

		// indexing from zero
		int j = 0;

		// starting from the bottom, compute the least upper bound for each row
		for (int i = p - 1; i >= 0; i--) {
			// search this row until we find a value less than u
			while (j < q && h.apply(i, j) > u) {
				j++;
			}
			// the entry preceding the one we just found is greater than u
			P[i] = j - 1;
		}

		return P;
	}

	private static int[] lessThanH(BiFunction<Integer, Integer, Double> h, int p, int q, double u) {

		// vector of size p
		final int[] Q = new int[p];

		// last possible row index
		int j = q - 1;

		// starting from the top, compute the greatest lower bound for each row
		for (int i = 0; i < p; i++) {
			// search this row until we find a value greater than u
			while (j >= 0 && h.apply(i, j) < u) {
				j--;
			}
			// the entry following the one we just found is less than u
			Q[i] = j + 1;
		}

		return Q;
	}

	/**
	 * Calculates the median for a list of values (<code>Number</code> objects). If <code>copyAndSort</code> is <code>false</code>, the list is assumed to be
	 * presorted in ascending order by value.
	 * 
	 * @param values
	 *            the values (<code>null</code> permitted).
	 * @param copyAndSort
	 *            a flag that controls whether the list of values is copied and sorted.
	 * 
	 * @return The median.
	 */
	public static double calculateMedian(TDoubleList values, boolean copyAndSort) {

		double result = Double.NaN;
		if (values != null) {
			if (copyAndSort) {
				final TDoubleList copy = new TDoubleArrayList(values);
				copy.sort();
				values = copy;
			}
			result = calcMedian(values, 0, values.size() - 1);
		}
		return result;
	}

	private static double calcMedian(TDoubleList values, int start, int end) {

		double result = Double.NaN;
		final int count = end - start + 1;
		if (count > 0) {
			if (count % 2 == 1) {
				if (count > 1) {
					result = values.get(start + ((count - 1) / 2));
				} else {
					result = values.get(start);
				}
			} else {
				final double value1 = values.get(start + ((count / 2) - 1));
				final double value2 = values.get(start + (count / 2));
				result = (value1 + value2) / 2.0;
			}
		}
		return result;
	}

	// Main call is WeightedMedian(a, 1, n) // a[1..n], p, r
	// Returns lower median
	private static double calcWeightedMedian(TDoubleList values, TDoubleList weights, int start, int end) {

		// base case for single element
		if (end == start) {
			return values.get(end);
		}
		// base case for two elements
		// make sure we return the average, in case the two candidates have equal weight
		if (end - start == 1) {
			if (Precision.equals(weights.get(start), weights.get(end))) {
				return (values.get(start) + values.get(end)) / 2.0;
			} else if (weights.get(start) > weights.get(end)) {
				return values.get(start);
			} else {
				return values.get(end);
			}
		}
		// partition around pivot r
		final int q = partition(a, p, r);
		final double wl = weights.subList(start, q).sum();
		final double wg = weights.subList(q + 1, end + 1).sum();
		// if partitions are balanced then we are done
		if (wl < 0.5 && wg < 0.5) {
			return values.get(q);
		} else {
			// increase pivot weight by the amount of partition we eliminate
			if (wl > wg) {
				weights.set(q, weights.get(q) + wg);
				// recurse on pivot inclusively
				return calcWeightedMedian(values, weights, start, q);
			} else {
				weights.set(q, weights.get(q) + wl);
				return calcWeightedMedian(values, weights, q, end);
			}
		}
	}

	public static double calculateMAD(TDoubleList values, double median) {

		double result = Double.NaN;
		if (values != null) {
			final TDoubleList diffs = new TDoubleArrayList(values.size());
			for (int i = 0; i < values.size(); i++) {
				diffs.add(Math.abs(values.get(i) - median));
			}
			diffs.sort();
			final int count = diffs.size();
			if (count > 0) {
				if (count % 2 == 1) {
					if (count > 1) {
						result = diffs.get((count - 1) / 2);
					} else {
						result = diffs.get(0);
					}
				} else {
					final double value1 = diffs.get(count / 2 - 1);
					final double value2 = diffs.get(count / 2);
					result = (value1 + value2) / 2.0;
				}
			}
		}
		return result;
	}

	public static double calculateMAD(TDoubleList values, boolean copyAndSort) {

		double result = Double.NaN;
		if (values != null) {
			final double median = calculateMedian(values, copyAndSort);
			result = calculateMAD(values, median);
		}
		return result;
	}

	/**
	 * Returns the standard deviation of a set of numbers.
	 * 
	 * @param data
	 *            the data (<code>null</code> or zero length array not permitted).
	 * 
	 * @return The standard deviation of a set of numbers.
	 */
	public static double getVariance(TDoubleList data, double mean) {

		double result = Double.NaN;
		if (data != null && data.size() > 0 && !Double.isNaN(mean)) {
			double sum = 0.0;
			for (int i = 0; i < data.size(); i++) {
				final double n = data.get(i);
				final double diff = n - mean;
				sum = sum + (diff * diff);
			}
			result = sum / (data.size());
		}
		return result;
	}

	public static double getVariance(TDoubleList data) {

		final double avg = calculateMean(data);
		return getVariance(data, avg);
	}

	public static double getStdDev(TDoubleList data) {

		final double avg = calculateMean(data);
		return getStdDev(data, avg);
	}

	public static double getStdDev(TDoubleList data, double mean) {

		return Math.sqrt(getVariance(data, mean));
	}

}
