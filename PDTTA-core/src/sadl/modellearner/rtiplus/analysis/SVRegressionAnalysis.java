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

import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_parameter;
import libsvm.svm_problem;

public class SVRegressionAnalysis extends DistributionAnalysis {

	public SVRegressionAnalysis(DistributionAnalysis fewElementsAnalysis, int fewElementsLimit) {
		super(fewElementsAnalysis, fewElementsLimit);
		// TODO Auto-generated constructor stub
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

		final svm_problem prob = new svm_problem();
		prob.y = new double[frequencies.size()];
		prob.l = values.size();
		prob.x = new svm_node[values.size()][1];

		for (int i = 0; i < values.size(); i++){
			final svm_node node = new svm_node();
			node.index = 1;
			node.value = normalize(values.get(i), valMin, valMax);
			prob.x[i][0] = node;
			prob.y[i] = normalize(frequencies.get(i), freqMin, freqMax);
		}

		final svm_parameter param = new svm_parameter();
		param.probability = 1;
		param.gamma = 0.5;
		param.nu = 0.5;
		param.C = 1000;
		param.svm_type = svm_parameter.EPSILON_SVR;
		param.kernel_type = svm_parameter.RBF;
		param.cache_size = 20000;
		param.eps = 0.001;

		final svm_model model = svm.svm_train(prob, param);

		for (int i = begin; i <= end; i++) {
			final svm_node node = new svm_node();
			node.index = 1;
			node.value = normalize(i, valMin, valMax);
			final double x = svm.svm_predict(model, new svm_node[] { node });
			System.out.println(i + "\t-> " + x);
		}



		// TODO Auto-generated method stub
		return new TIntArrayList();
	}

	public static void main(String[] args) {

		final SVRegressionAnalysis a = new SVRegressionAnalysis(null, -1);

		final TIntList v = new TIntArrayList(new int[] { 1, 5, 6, 7, 8, 20, 21, 22, 23, 24, 50, 99 });
		final TIntList f = new TIntArrayList(new int[] { 1, 1, 10, 1, 1, 10, 1, 1, 10, 10, 10, 99 });

		// final List<TIntList> y = a.computeClusters(v, f);
		// System.out.println("Clusters: " + y.toString());

		final TIntList x = a.performAnalysis(v, f, 0, 100);
		System.out.println("Splits: " + x.toString());

	}

}
