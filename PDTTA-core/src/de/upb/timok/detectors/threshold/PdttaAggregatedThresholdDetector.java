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
package de.upb.timok.detectors.threshold;

import gnu.trove.list.TDoubleList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import de.upb.timok.constants.ProbabilityAggregationMethod;
import de.upb.timok.detectors.PdttaDetector;

public class PdttaAggregatedThresholdDetector extends PdttaDetector {
	private static Logger logger = LoggerFactory.getLogger(PdttaAggregatedThresholdDetector.class);

	double aggregatedEventThreshold;
	double aggregatedTimeThreshold;

	public PdttaAggregatedThresholdDetector(ProbabilityAggregationMethod aggType, double aggregatedEventThreshold, double aggregatedTimeThreshold) {
		super(aggType);
		this.aggregatedEventThreshold = aggregatedEventThreshold;
		this.aggregatedTimeThreshold = aggregatedTimeThreshold;
	}

	@Override
	protected boolean decide(TDoubleList eventLikelihoods, TDoubleList timeLikelihoods) {
		final double normalizedEventThreshold = aggregatedEventThreshold;
		final double normalizedTimeThreshold = aggregatedTimeThreshold;
		// TODO: Was passiert, wenn aggregate NAN zurück gibt?!
		// Sollte eigentlich nur -infty zurückgeben
		final double aggregatedEventScore = aggregate(eventLikelihoods, aggType);
		final double aggregatedTimeScore = aggregate(timeLikelihoods, aggType);
		logger.info("aggEventScore={}\taggTimeScore={}", aggregatedEventScore, aggregatedTimeScore);

		if (aggregatedEventScore <= normalizedEventThreshold || aggregatedTimeScore <= normalizedTimeThreshold) {
			// one of the scores is so low, that this indicates an anomaly
			return true;
		} else {
			return false;
		}
	}

}
