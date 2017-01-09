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

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

public class PolynomialRegressionAnalysis extends DistributionAnalysis {

	private final int degree;

	public PolynomialRegressionAnalysis(int degree) {
		this(degree, null, -1);
	}

	public PolynomialRegressionAnalysis(int degree, DistributionAnalysis fewElementsAnalysis, int fewElementsLimit) {
		super(fewElementsAnalysis, fewElementsLimit);

		this.degree = degree;
	}

	@Override
	TIntList analyzeDistribution(TIntList values, TIntList frequencies, int begin, int end) {

		int valMin = Integer.MAX_VALUE;
		int valMax = Integer.MIN_VALUE;
		final TIntIterator it = values.iterator();
		while (it.hasNext()) {
			final int v = it.next();
			if (v < valMin) {
				valMin = v;
			}
			if (v > valMax) {
				valMax = v;
			}
		}

		int freqMin = Integer.MAX_VALUE;
		int freqMax = Integer.MIN_VALUE;
		final TIntIterator it2 = values.iterator();
		while (it2.hasNext()) {
			final int v = it2.next();
			if (v < freqMin) {
				freqMin = v;
			}
			if (v > freqMax) {
				freqMax = v;
			}
		}

		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int i = 0; i < values.size(); i++) {
			obs.add(normalize(values.get(i), valMin, valMax), normalize(frequencies.get(i), freqMin, freqMax));
		}
		final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);
		final double[] coeff = fitter.fit(obs.toList());

		final LaguerreSolver laguerreSolver = new LaguerreSolver();
		final Complex[] roots = laguerreSolver.solveAllComplex(coeff, begin - 1);

		final int min = valMin;
		final int max = valMax;
		final int[] orgRoots = Arrays.stream(roots).mapToInt(r -> denormalize(r.getReal(), min, max)).distinct().sorted().toArray();
		return new TIntArrayList(orgRoots);
	}

	public static void main(String[] args) {

		final PolynomialRegressionAnalysis a = new PolynomialRegressionAnalysis(10);

		final TIntList v = new TIntArrayList(new int[] { 1, 5, 6, 7, 8, 20, 21, 22, 23, 24, 50, 99 });
		final TIntList f = new TIntArrayList(new int[] { 1, 1, 10, 1, 1, 10, 1, 1, 10, 10, 10, 99 });

		// final List<TIntList> y = a.computeClusters(v, f);
		// System.out.println("Clusters: " + y.toString());

		final TIntList x = a.performAnalysis(v, f, 0, 100);
		System.out.println("Splits: " + x.toString());

	}

}
