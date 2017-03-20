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
package sadl.modellearner.rtiplus.analysis;

import gnu.trove.list.TDoubleList;
import sadl.modellearner.rtiplus.StatisticsUtil;

/**
 * Adjusted IQR Analysis for skewed distributions</br>
 * </br>
 * Ref.: M. Hubert, E. Vandervieren, An Adjusted Boxplot for Skewed Distributions, Computational Statistics & Data Analysis, Volume 52, Issue 12, 15 August
 * 2008, Pages 5186-5201</br>
 * <a href="https://wis.kuleuven.be/stat/robust/papers/2008/adjboxplot-revision.pdf">PDF</a>
 * 
 * @author Fabian Witter
 *
 */
public class AdjustedIQROutlierAnalysis extends OutlierDistanceAnalysis {

	private final boolean onlyFarOuts;
	private final boolean presumeOnlyFewOutliers;

	public AdjustedIQROutlierAnalysis(double strength, boolean onlyFarOuts, boolean presumeOnlyFewOutliers) {
		this(strength, onlyFarOuts, presumeOnlyFewOutliers, null, -1);
	}

	public AdjustedIQROutlierAnalysis(double strength, boolean onlyFarOuts, boolean presumeOnlyFewOutliers, DistributionAnalysis fewElementsAnalysis,
			int fewElementsLimit) {
		super(strength, fewElementsAnalysis, fewElementsLimit);
		this.onlyFarOuts = onlyFarOuts;
		this.presumeOnlyFewOutliers = presumeOnlyFewOutliers;
	}

	@Override
	int getOutlierDistance(TDoubleList distValues) {

		distValues.sort();
		final double q1 = StatisticsUtil.calculateQ1(distValues, false);
		final double q3 = StatisticsUtil.calculateQ3(distValues, false);
		final double median = StatisticsUtil.calculateMedian(distValues, false);
		distValues.reverse();
		final double medcouple = StatisticsUtil.calculateMedcouple(distValues, median, false);
		double coeff = onlyFarOuts ? 3.0 : 1.5;
		if (presumeOnlyFewOutliers) {
			// Do not care about direction of skewness
			coeff *= Math.exp(3.5 * medcouple);
		} else if (medcouple > 0) {
			// Skewed to the right
			coeff *= Math.exp(3 * medcouple);
		} else {
			// Skewed to the left
			coeff *= Math.exp(4 * medcouple);
		}
		return (int) Math.ceil(((q3 + (q3 - q1) * coeff) / 2.0));
	}

}
