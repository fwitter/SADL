package sadl.modellearner.rtiplus.analysis;

import org.apache.commons.math3.util.Precision;

// The software is made available under the MIT license:

//Copyright (c) 2015 Andre Rauh
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//THE SOFTWARE.

// https://github.com/rauhs/fast-weighted-median/blob/master/src/wmedianf_impl.cpp

public class MedianOutlierAnalysis {

	// Factor by which the left side is multiplied to yield the right side when
	// using the bisection method to get the optimal k for a given M
	static final double BISECTION_DELTA = 1.1;
	// Tolerance at which to stop when using the bisection method
	static final double BISECTION_TOL = 1e-2;
	// If alph*M is smaller than this value the algorithm will use poisson
	// distribution to model the cost function. Otherwise normal distribution will hold.
	static final int NORMAL_POISSON_LIMIT = 10;

	// Finds the median of the elements x[0], x[dist], x[2dist] ... x[(N-1)*dist]
	static double medianf_equidist(double[] a, double[] w, int l, int r, int k, int dist) {
		int n, i;
		int j;
		int s, sd;
		int ll, rr;
		double z;
		double t;

		if (dist <= 0) {
			dist = 1;
		}

		while (r > l) {
			if ((r - l) > 600) {
				n = r - l + 1;
				i = k - l + 1;
				z = Math.log(n);
				s = (int) (0.5 * Math.exp(2 * z / 3));
				sd = (int) (0.5 * Math.sqrt(z * s * (n - s) / n));
				if (i - n / 2.0 < 0) {
					sd = -sd;
				}
				ll = (int) Math.max(l, k - (double) i * s / n + sd);
				rr = (int) Math.min(r, k + (double) (n - i) * s / n + sd);
				medianf_equidist(a, w, ll, rr, k, dist);
			}
			t = a[k * dist]; // Pivot
			i = l;
			j = r;
			swap(a, l * dist, k * dist);
			swap(w, l * dist, k * dist);
			if (a[r * dist] > t) {
				swap(a, r * dist, l * dist);
				swap(w, r * dist, l * dist);
			}
			while (i < j) {
				swap(a, i * dist, j * dist);
				swap(w, i * dist, j * dist);
				i++;
				j--;
				while (a[i * dist] < t) {
					i++;
				}
				while (t < a[j * dist]) {
					j--;
				}
			}
			if (a[l * dist] == t) {
				swap(a, l * dist, j * dist);
				swap(w, l * dist, j * dist);
			} else {
				j++;
				swap(a, j * dist, r * dist);
				swap(w, j * dist, r * dist);
			}
			if (j <= k) {
				l = j + 1;
			}
			if (k <= j) {
				r = j - 1;
			}
		}
		return a[k * dist];
	}

	static double wmedianf(double[] x, double[] w, double initialPivot, double W0) {

		final int K = 32; // If problem size is smaller than this abort and find median using quickselect
		int M0; // Number of samples to get the first pivot
		int left = 0; // Always points to the first element of (reduced) problem
		int right = x.length - 1; // Always points to the last element of (reduced) problem
		double pivot = x[0];
		double xmax = x[0];
		double xmin = x[0];
		final double Wl = 0.0; // Sum of weight discarded to the left (0...left-1)
		final double Wr = 0.0; // Sum of weight discarded to the right (right-1...END)

		// Stop recursion
		if (x.length < 256) {
			return wquickSelect(x, w, x.length, W0, Wl, Wr);
		}

		// Calculate the number of samples to get the first pivot
		M0 = getM(x.length);

		// If the passed pivot is null pointer then compute the first pivot
		if (Double.isNaN(initialPivot)) {
			// Need to choose the first pivot
			pivot = medianf_equidist(x, w, 0, M0 - 1, (M0 - 1) / 2, x.length / M0);

			// Partitioning will fail if pivot happend to be NaN
			if (Double.isNaN(pivot)) {
				return wquickSelect(x, w, x.length, W0, Wl, Wr);
			}
		} else {
			pivot = initialPivot;
		}

		// partition will return the index of the first element >(=) than the pivot
		int lrb = 0;// Indicates to count both and compute W0
		final int idx0 = wpartition(x, w, pivot, x.length, lrb, W0, Wl, Wr);

		if (lrb > 0) {
			// Pivot was larger than the median
			right = idx0;
			xmax = pivot;
		} else {
			// Pivot was smaller than the median
			left = idx0;
			xmin = pivot;
		}
		checkConsistency(x, w, W0, Wl, Wr, left, right, x.length);

		// Note: In case the pivot was the WM we just continue until the end where
		// the standard Quickselect solves the remaining problem.

		int N1 = right - left + 1; // Number of samples left

		// While the problem isn't bounded above and below:
		while ((left == 0) || (right == x.length - 1)) {
			double alph;
			int M1;

			if (N1 < K) {
				return wquickSelect(x[left], w[left], N1, W0, Wl, Wr);
			}

			// This can occur if negative weights are given.
			if ((Wr < 0) || (Wl < 0)) {
				throw new IllegalArgumentException("No negative weigths allowed");
			}

			if (left == 0) {
				alph = W0 / (2 * W0 - Wr);
				M1 = getScaledM(M0, (W0 - 0.5 * Wr) / W0);
			} else {
				alph = 1 - W0 / (2 * W0 - Wl);
				M1 = getScaledM(M0, (W0 - 0.5 * Wl) / W0);
			}
			if (alph <= 0 || alph >= 1) {
				// We found the median
				return pivot;
			}
			// Let M1 never be less than 5:
			M1 = Math.max(M1, 5);
			if (M1 > N1) {
				M1 = N1;
			}

			int k = getK(M1, alph);

			// Never choose min/max as pivot
			k = Math.max(3, k); // Never select the min
			k = Math.min(M1 - 2, k); // Never select the max

			pivot = medianf_equidist(x[left], w[left], 0, M1 - 1, k - 1, N1 / M1);
			// Partitioning will fail if pivot happend to be NaN
			if (Double.isNaN(pivot)) {
				return wmSortNanAware(x[left], w[left], N1, W0 - Wl);
			}

			// partition will return the index of the first element >(=) to the pivot
			lrb = lrb; // lrb stays the same so to count the hopefully smaller side
			checkConsistency(x, w, W0, Wl, Wr, left, right, x.length);
			int idx1 = wpartition(x[left], w[left], pivot, N1, lrb, W0, Wl, Wr);

			idx1 += left; // partition() will count from x[left]
			if (lrb > 0) {
				// Pivot was larger than the median
				right = idx1;
				xmax = pivot;
			} else {
				// Pivot was smaller than the median
				left = idx1;
				xmin = pivot;
			}
			checkConsistency(x, w, W0, Wl, Wr, left, right, x.length);
			N1 = right - left + 1; // Update samples left
		}
		int N2plus = N1;
		int elements_removed = 0;

		// Start reducing the set using the approximation of the median
		for (;;) {
			if (N2plus <= K) {
				break;
			}
			// If any of the prev pivots were inf/NaN then we cannot take a lin
			// combination of them and have to fall back.
			if (Double.isInfinite(xmin) || Double.isInfinite(xmax)) {
				break;
			}

			// 0 < a < 1
			final double a = (W0 - Wl) / (2 * W0 - Wl - Wr);
			final double c = approx_fct(a);
			// This pivot is problably not in the set.
			pivot = c * xmax + (1 - c) * xmin;

			if (c < 0.5) {
				lrb = -1;
			} else {
				lrb = +1;
			}
			int idx2plus = wpartition(x[left], w[left], pivot, N2plus, lrb, W0, Wl, Wr);

			idx2plus += left; // wpartition() will count from x[left]
			if (lrb > 0) {
				// Pivot was larger than the median
				right = idx2plus;
				xmax = pivot;
			} else {
				// Pivot was smaller than the median
				left = idx2plus;
				xmin = pivot;
			}
			checkConsistency(x, w, W0, Wl, Wr, left, right, x.length);
			elements_removed = N2plus - (right - left + 1);
			N2plus = right - left + 1;
			if (elements_removed <= 2) {
				// We're stuck -> fallback to old method
				break;
			}
		}
		// variable median is one-based index
		checkConsistency(x, w, W0, Wl, Wr, left, right, x.length);
		return wquickSelect(x[left], w[left], N2plus, W0, Wl, Wr);
	}

	// Note: pivot is not necessarily in the set
	// lrb:
	// - (-1): count from left
	// - (+1): count from right
	// - all other: count both
	// WO is only written if lrb is not +-1.
	// lrb is also a return value and is set to either -1 or +1:
	// o -1 indicates that the pivot was less than the WM
	// o +1 indicates that the pivot was less than the WM
	// Wl and Wr will be set to the sum of the weigts of the elements removed.
	static int wpartition(double[] x, double[] w, double pivot, int N, int lrb, double W0, double Wl, double Wr) {
		// summed up concomitant weights which are smaller/larger than the pivot
		double wleftSum = 0.0;
		double wrightSum = 0.0;

		int left = 0;
		int right = N - 1;

		if (lrb == -1) {
			// W0 must be provided if lrb is -1
			if (W0 <= 0.0) {
				throw new IllegalArgumentException("No W0 provided!");
			}

			// Sum up the weights left of the pivot
			while (left <= right) {
				while (x[left] < pivot) {
					wleftSum += w[left];
					left++;
				}
				while (x[right] > pivot) {
					right--;
				}

				if (left <= right) {
					swap(x, left, right);
					swap(w, left, right);
					wleftSum += w[left];
					left++;
					right--;
				}
			}
			// In case there were elements equal to the pivot, we need to adjust the sum of weights
			while (left > right + 1) {
				left--;
				wleftSum -= w[left];
			}
			if (wleftSum + Wl < W0) {
				Wl += wleftSum;
				lrb = -1;
				return left;
			} else {
				Wr = 2 * W0 - wleftSum - Wl;
				lrb = +1;
				return right;
			}
		} else if (lrb == 1) {
			if (W0 <= 0.0) {
				throw new IllegalArgumentException("No W0 provided!");
			}

			// Sum up the weights right of the pivot
			while (left <= right) {
				while (x[left] < pivot) {
					left++;
				}
				while (x[right] > pivot) {
					wrightSum += w[right];
					right--;
				}

				if (left <= right) {
					swap(x, left, right);
					swap(w, left, right);
					wrightSum += w[right];
					++left;
					--right;
				}
			}
			// In case there were elements equal to the pivot, we need to adjust the sum of weights
			while (left > right + 1) {
				left--;
				wleftSum -= w[left];
			}
			if (wrightSum + Wr < W0) {
				Wr += wrightSum;
				lrb = +1;
				return right;
			} else {
				Wl = 2 * W0 - wrightSum - Wr;
				lrb = -1;
				return left;
			}
		} else {
			// Sum up the weights right of the pivot
			while (left <= right) {
				while (x[left] < pivot) {
					wleftSum += w[left];
					left++;
				}
				while (x[right] > pivot) {
					wrightSum += w[right];
					right--;
				}
				if (left <= right) {
					swap(x, left, right);
					swap(w, left, right);
					wleftSum += w[left];
					wrightSum += w[right];
					++left;
					--right;
				}
			}
			// In case there were elements equal to the pivot, we need to adjust the sum of weights
			while (left > right + 1) {
				left--;
				wleftSum -= w[left];
			}
			W0 = 0.5 * (wleftSum + wrightSum);
			if (wleftSum < W0) {
				Wl = wleftSum;
				lrb = -1;
				return left;
			} else {
				Wr = wrightSum;
				lrb = +1;
				return right;
			}
		}
	}

	/*
	 * This Quickselect routine is based on the algorithm described in
	 * "Numerical recipes in C", Second Edition,
	 * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	 * This code by Nicolas Devillard - 1998. Public domain.
	 */
	static double wquickSelect(double[] x, double[] w, int N, double W0, double Wl, double Wr) {

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
			if (high <= low) /* One element only */
			{
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

	// Returns the optimal M to a given N
	// Where N is the size of the array and M is the size of the subset
	static int getM(int N) {
		int M = (int) (Math.pow(N, 0.618764) - Math.sqrt(N) + 5.0);
		M |= 0x01; // Make M odd
		M = Math.max(3, M); // Make M at least 3
		return M;
	}

	// Returns the optimal new M to a given old M and a factor
	static int getScaledM(int M, double scale) {
		return (int) (M * Math.pow(scale, 0.618764));
	}

	// Returns the optimal OS (k) for a given alpha & M
	static int getK(int M_in, double alph_in) {
		double alph = alph_in;
		final double Md = M_in;

		// Make alph always 0..1/2
		if (alph > 0.5) {
			alph = 1 - alph_in;
		}

		// Fitting was done with data a in 1/2..1
		final double a = 1 - alph;
		// This is the function from Eureqa (MSE ~ 500) (fits from M = 50..12000)
		double xold = 0.998202 * a * M_in + 11.421 * M_in / Math.pow(92.9553 * M_in, a) - 1.86877 * Math.pow(a * M_in - a * a * M_in, 0.545273);
		xold = M_in - xold; // Invert since a is in 0..1/2

		// Ad hoc adjustment for a starting point left of the optimal point
		xold = 0.5 * alph * M_in + 0.5 * xold;

		if (M_in * alph < NORMAL_POISSON_LIMIT) {
			// For poisson use lower starting point to ensure convergence
			xold = 0.6 * alph * M_in + 0.4 * xold;
			// Make starting point never < 1 -> This helps extremly small values of alph
			if (xold < 1) {
				xold = 1;
			}
		}

		// The final value will go here
		double xnew = 0.0;

		if (M_in < 12000) {
			// Fitting holds and Newton method will converge
			xnew = xold - fdivfp(xold, alph, M_in);
			xold = xnew;
			xnew = xold - fdivfp(xold, alph, M_in);
		} else {
			// M > 12000 -> Use bisection method
			double xleft = alph * M_in;
			double xright = BISECTION_DELTA * xleft;
			double yleft = cost_fct(xleft, alph, M_in);
			double yright = cost_fct(xright, alph, M_in);
			// First make sure that yleft & yright have different signs
			while (yright <= 0) {
				xleft = xright;
				yleft = yright; // Make the left point the right point
				xright = BISECTION_DELTA * xleft;
				yright = cost_fct(xright, alph, M_in);
			}
			// Now we ensured that yright is positive.
			if (xleft <= 0) {
				// This should never happen
				throw new IllegalArgumentException("xleft <= 0");
			}

			int j = 0;
			while ((xright - xleft) / M_in > BISECTION_TOL) {
				// Take the midpoint
				final double xmid = (xleft + xright) / 2;
				final double ymid = cost_fct(xmid, alph, M_in);
				if (ymid <= 0) {
					// It's on the left side of the cost function
					xleft = xmid;
					yleft = ymid;
				} else {
					xright = xmid;
					yright = ymid;
				}
				if (++j > 8) {
					break; // Do at most 8 iterations
				}
			}
			// The final output is going to be the midpoint
			xnew = (xleft + xright) / 2;
		}

		xnew += 1.0;

		// Check if first pivot was left/right of median.
		if (alph_in > 0.5) {
			xnew = M_in - xnew + 1.0;
		}

		return (int) xnew;
	}

	static double cost_fct(double x, double alph, int M_in) {
		if (M_in * alph >= NORMAL_POISSON_LIMIT) {
			// Use normal approx
			return (2 * x - M_in + 1) * normpdf(x, alph * M_in, alph * (1 - alph) * M_in) + 1 - 2 * betai(x + 1, M_in - x, alph);
		} else {
			// Use poisson approx
			return (2 * x - M_in + 1) * poisspdf(x, alph * M_in) + 1;
		}
	}

	static double fdivfp(double x, double alph, int M_in) {
		final double mu = M_in * alph;
		double num, den;

		if (mu >= NORMAL_POISSON_LIMIT) {
			// Use normal approx
			// fdivfpN = @(k) (2*k-M+1 + (1-2*betainc(alph(Nidx),k+1,M-k))./normpdf(k,alph(Nidx)*M,sig(Nidx))) ...
			// ./ (4+(2*k-M+1).*(alph(Nidx)*M-k)./sig(Nidx).^2);
			final double sigma = Math.sqrt(mu * (1 - alph));
			num = 2 * x - M_in + 1 + (1 - 2 * betai(x + 1, M_in - x, alph)) / normpdf(x, mu, sigma);
			den = 4 + (2 * x - M_in + 1) * (mu - x) / (sigma * sigma);
		} else {
			// Use poisson approx
			// fdivfpP = @(k) (k - (M-1)/2 + 1/2./poisspdf(k,M*alph(Pidx))) ...
			// ./ (1 + (k - (M-1)/2).*(log(M*alph(Pidx)./k) - 1./(2*k) + 1/12./k.^2));
			final double mu2 = M_in * alph;
			final double xn = x - 0.5 * (M_in - 1);
			num = xn + 0.5 / poisspdf(x, mu2);
			den = 1 + xn * Math.log(mu2 / x) - 1 / (2 * x);
		}
		return num / den;
	}

	static double approx_fct_1(double x) {
		return x - x * x + x * x * x + 1.4 * x * x * x * x;
	}

	/**
	 * Interpolation function used for calculating the pivot when the method of
	 * linear combination of max & min is used.
	 */
	static double approx_fct(double a) {
		if (a < .5) {
			return .5 - approx_fct_1(.5 - a);
		} else {
			return .5 + approx_fct_1(-.5 + a);
		}
	}

	static void swap(double[] x, int a, int b) {
		final double temp = x[a];
		x[a] = x[b];
		x[b] = temp;
	}

	/**
	 * Solving the WM problem by sorting. If W0 is unknown it can be passed a negative
	 * value and the function will compute it.
	 */
	static double wmSort(double[] x, double[] w, int N, double W0) {
		ShellSortPair(x, w, N);

		if (W0 < 0.0) {
			W0 = 0.0;
			for (int i = 0; i < N; i++) {
				W0 += w[i];
			}
			W0 *= 0.5;
		}

		double wSum = 0.0;
		int i;

		for (i = 0; i < N; i++) {
			wSum += w[i];
			if (wSum >= W0) {
				break;
			}
		}
		return x[i];
	}

	/**
	 * A shell sorting routine which works on x and shuffles the pair (x,w) accordingly
	 */
	static void ShellSortPair(double[] x, double[] w, int size) {

		final int hmax = size / 9;
		int h;
		for (h = 1; h <= hmax; h = 3 * h + 1) {
			;
		}
		for (; h > 0; h /= 3) {
			for (int i = h; i < size; ++i) {
				final double v = x[i];
				final double vw = w[i];
				int j = i;
				while (j >= h && v < x[j - h]) {
					x[j] = x[j - h];
					w[j] = w[j - h];
					j -= h;
				}
				x[j] = v;
				w[j] = vw;
			}
		}
	}

	/**
	 * Floyd and Rivest SELECT algorithm.
	 * // left is the left index for the interval
	 * // right is the right index for the interval
	 * // k is the desired index value, where array[k] is the (k+1)th smallest element when left = 0
	 */
	public static double select_floyd_rivest(double[] a, int l, int r, int k) {

		int n, i, j, s, sd, ll, rr;
		double z, t;

		while (r > l) {
			// use select recursively to sample a smaller set of size s
			// the arbitrary constants 600 and 0.5 are used in the original
			// version to minimize execution time
			if ((r - l) > 600) {
				n = r - l + 1;
				i = k - l + 1;
				z = Math.log(n);
				s = (int) (0.5 * Math.exp(2 * z / 3));
				sd = (int) (0.5 * Math.sqrt(z * s * (n - s) / n));
				if (i - n / 2.0 < 0) {
					sd = -sd;
				}
				ll = (int) Math.max(l, k - (double) i * s / n + sd);
				rr = (int) Math.min(r, k + (double) (n - i) * s / n + sd);
				select_floyd_rivest(a, ll, rr, k);
			}
			// partition the elements between left and right around t
			t = a[k];
			i = l;
			j = r;
			swap(a, l, k);
			if (a[r] > t) {
				swap(a, r, l);
			}
			while (i < j) {
				swap(a, i, j);
				i++;
				j--;
				while (a[i] < t) {
					i++;
				}
				while (t < a[j]) {
					j--;
				}
			}
			if (a[l] == t) {
				swap(a, l, j);
			} else {
				j++;
				swap(a, j, r);
			}
			// adjust left and right towards the boundaries of the subset
			// containing the (k - left + 1)th smallest element
			if (j <= k) {
				l = j + 1;
			}
			if (k <= j) {
				r = j - 1;
			}
		}
		return a[k];
	}

	/**
	 * Floyd and Rivest SELECT algorithm.
	 * Takes an additional parameter delta which is an integer specifying the inter sample
	 * distance
	 */
	double select_floyd_rivest_delta(double[] a, int l, int r, int k, int delta) {
		int n, i, j, s, sd, ll, rr;
		double z, t;

		if (((k - l) % delta != 0) || ((r - l) % delta != 0)) {
			return -1.0; // assume k is l + some multiple of delta
		}

		while (r > l) {
			if ((r - l) > 600 * delta) {
				n = (r - l + 1) / delta;
				i = (k - l + 1) / delta;
				z = Math.log(n);
				s = (int) (0.5 * Math.exp(2 * z / 3));
				sd = (int) (0.5 * Math.sqrt(z * s * (n - s) / n));
				if (i - n / 2.0 < 0) {
					sd = -sd;
				}
				ll = Math.max(l, k - ((int) ((double) i * s / n + sd)) * delta);
				rr = Math.min(r, k + ((int) ((double) (n - i) * s / n + sd)) * delta);
				select_floyd_rivest(a, ll, rr, k);
			}
			t = a[k];
			i = l;
			j = r;
			swap(a, l, k);
			if (a[r] > t) {
				swap(a, r, l);
			}
			while (i < j) {
				swap(a, i, j);
				i += delta;
				j -= delta;
				while (a[i] < t) {
					i += delta;
				}
				while (t < a[j]) {
					j -= delta;
				}
			}
			if (a[l] == t) {
				swap(a, l, j);
			} else {
				j += delta;
				swap(a, j, r);
			}
			if (j <= k) {
				l = j + delta;
			}
			if (k <= j) {
				r = j - delta;
			}
		}
		return a[k];
	}

	/**
	 * A debug routine to check if the values Wl, W0, Wr etc. still are in agreement
	 * with each other
	 */
	static void checkConsistency(double[] x, double[] w, double W0_ref, double Wl_ref, double Wr_ref, int left, int right, int N) {
		// #if _DEBUG TODO
		double Wl = 0.0;
		double Wr = 0.0;
		double W0 = 0.0;

		int i = 0;
		for (; i < left; i++) {
			Wl += w[i];
		}
		for (; i <= right; i++) {
			W0 += w[i];
		}
		for (; i < N; i++) {
			Wr += w[i];
		}
		W0 += Wr + Wl;
		W0 *= 0.5;

		if (!Precision.equals(Wl, Wl_ref)) {
			throw new RuntimeException("Wl not correct");
		}
		if (!Precision.equals(Wr, Wr_ref)) {
			throw new RuntimeException("Wr not correct");
		}
		if (!Precision.equals(W0, W0_ref)) {
			throw new RuntimeException("W0 not correct");
		}
	}

}