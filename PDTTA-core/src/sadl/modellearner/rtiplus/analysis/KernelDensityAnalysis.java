/**
 * This file is part of SADL, a library for learning all sorts of (timed) automata and performing sequence-based anomaly detection.
 * Copyright (C) 2013-2017  the original author or authors.
 *
 * SADL is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * SADL is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with SADL.  If not, see <http://www.gnu.org/licenses/>.
 */
package sadl.modellearner.rtiplus.analysis;

import java.util.Arrays;

import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.linked.TDoubleLinkedList;
import gnu.trove.list.linked.TIntLinkedList;
import jsat.distributions.empirical.KernelDensityEstimator;
import jsat.distributions.empirical.kernelfunc.GaussKF;
import jsat.distributions.empirical.kernelfunc.KernelFunction;
import jsat.linear.DenseVector;
import jsat.linear.Vec;

public class KernelDensityAnalysis extends DistributionAnalysis {

	private final KernelFunction kf;
	private final double bandwidth;

	public KernelDensityAnalysis(KernelFunction kernelFkt, double KdeBandwidth) {
		this(kernelFkt, KdeBandwidth, null, -1);
	}

	public KernelDensityAnalysis(KernelFunction kernelFkt, double KdeBandwidth, DistributionAnalysis fewElementsAnalysis, int fewElementsLimit) {
		super(fewElementsAnalysis, fewElementsLimit);

		this.kf = kernelFkt;
		this.bandwidth = KdeBandwidth;
	}

	@Override
	TIntList analyzeDistribution(TIntList values, TIntList frequencies, int begin, int end) {

		final TDoubleList transitionTimes = new TDoubleLinkedList();
		for (int i = 0; i < values.size(); i++) {
			final double[] vals = new double[frequencies.get(i)];
			Arrays.fill(vals, values.get(i));
			transitionTimes.add(vals);
		}

		final Vec v = new DenseVector(transitionTimes.toArray());

		final KernelDensityEstimator kde = new KernelDensityEstimator(v, kf, bandwidth);

		final TIntList splits = new TIntLinkedList();
		boolean gap = false;
		for (int i = begin - 1; i <= end; i++) {
			final double p = kde.pdf(i);
			if (p <= 0.0 && !gap) {
				splits.add(i);
				gap = true;
			} else if (p > 0.0 && gap) {
				splits.add(i);
				gap = false;
			}
		}

		return splits;
	}

	public static void main(String[] args) {

		final KernelDensityAnalysis a = new KernelDensityAnalysis(GaussKF.getInstance(), 0.1);

		final TIntList v = new TIntArrayList(new int[] { 1, 5, 6, 7, 8, 20, 21, 22, 23, 24, 50, 99 });
		final TIntList f = new TIntArrayList(new int[] { 1, 1, 10, 1, 1, 10, 1, 1, 10, 10, 10, 99 });

		// final List<TIntList> y = a.computeClusters(v, f);
		// System.out.println("Clusters: " + y.toString());

		final TIntList x = a.performAnalysis(v, f, 0, 100);
		System.out.println("Splits: " + x.toString());

	}

}
