/*******************************************************************************
 * This file is part of PDTTA, a library for learning Probabilistic deterministic timed-transition Automata.
 * Copyright (C) 2013-2015  Timo Klerx
 * 
 * PDTTA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * PDTTA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with PDTTA.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
package de.upb.timok.detectors;

import de.upb.timok.constants.ProbabilityAggregationMethod;
import de.upb.timok.constants.ScalingMethod;
import de.upb.timok.detectors.featureCreators.FeatureCreator;
import de.upb.timok.interfaces.TrainableDetector;
import de.upb.timok.oneclassclassifier.LibSvmClassifier;

@Deprecated
public class PdttaOneClassSvmDetector extends PdttaVectorDetector implements TrainableDetector {
	/**
	 * Use the PdttaVectorDetector constructor directly instead
	 * 
	 * @param aggType
	 * @param fc
	 * @param useProbability
	 * @param gamma
	 * @param nu
	 * @param costs
	 * @param kernelType
	 * @param eps
	 * @param degree
	 * @param scalingMethod
	 */
	public PdttaOneClassSvmDetector(ProbabilityAggregationMethod aggType, FeatureCreator fc, int useProbability, double gamma, double nu, double costs,
			int kernelType, double eps, int degree, ScalingMethod scalingMethod) {
		super(aggType, fc, new LibSvmClassifier(useProbability, gamma, nu, costs, kernelType, eps, degree, scalingMethod));
	}

}
