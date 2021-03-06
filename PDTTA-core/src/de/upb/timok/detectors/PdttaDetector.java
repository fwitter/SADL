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

import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import jsat.distributions.Distribution;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import de.upb.timok.constants.ProbabilityAggregationMethod;
import de.upb.timok.interfaces.AnomalyDetector;
import de.upb.timok.interfaces.Model;
import de.upb.timok.models.PDTTA;
import de.upb.timok.run.GenericSmacPipeline;
import de.upb.timok.structure.TimedSequence;
import de.upb.timok.structure.Transition;
public abstract class PdttaDetector implements AnomalyDetector {
	protected ProbabilityAggregationMethod aggType;
	private static Logger logger = LoggerFactory.getLogger(PdttaDetector.class);

	PDTTA model;
	@Override
	public boolean isAnomaly(Model newModel, TimedSequence s) {
		setModel(newModel);
		return isAnomaly(s);
	}

	@Override
	public boolean[] areAnomalies(Model newModel, List<TimedSequence> testSequences) {
		setModel(newModel);
		return areAnomalies(testSequences);

	}

	public PdttaDetector(ProbabilityAggregationMethod aggType) {
		super();
		this.aggType = aggType;
	}

	public PdttaDetector(ProbabilityAggregationMethod aggType, PDTTA model) {
		super();
		this.aggType = aggType;
		this.model = model;
	}

	/**
	 * returns two double values for every timed sequence. The first value is the event likelihood, the second the time likelihood
	 * 
	 * @param testTimedSequences
	 * @param transitionDistributions
	 * @return
	 */
	public List<double[]> computeAggregatedLikelihoods(List<TimedSequence> testTimedSequences) {
		final List<double[]> result = new ArrayList<>();
		for (final TimedSequence ts : testTimedSequences) {
			final double eventProb = computeAggregatedEventLikelihood(ts);
			final double timeProb = computeAggregatedTimeLikelihood(ts);
			result.add(new double[] { eventProb, timeProb });
		}
		return result;
	}

	/**
	 * returns two double values for every timed sequence. The first value is the event likelihood, the second the time likelihood
	 * 
	 * @param testTimedSequences
	 * @param transitionDistributions
	 * @return
	 */
	public double[] computeAggregatedLikelihood(TimedSequence testTimedSequence) {
		final double eventProb = computeAggregatedEventLikelihood(testTimedSequence);
		final double timeProb = computeAggregatedTimeLikelihood(testTimedSequence);
		return new double[] { eventProb, timeProb };
	}

	protected TDoubleList computeTimeLikelihoods(TimedSequence ts) {
		final TDoubleList list = new TDoubleArrayList();
		int currentState = 0;
		for (int i = 0; i < ts.getEvents().size(); i++) {
			final Transition t = model.getTransition(currentState, ts.getEvents().get(i));
			// DONE this is crap, isnt it? why not return an empty list or null iff there is no transition for the given sequence? or at least put a '0' in the
			// last slot.
			if (t == null) {
				list.add(0);
				return list;
			}
			final Distribution d = model.getTransitionDistributions().get(t.toZeroProbTransition());
			if (d == null) {
				// System.out.println("Found no time distribution for Transition "
				// + t);
				list.add(0);
			} else {
				list.add(d.pdf(ts.getTimeValues().get(i)));
			}
			currentState = t.getToState();
		}
		return list;
	}

	public double computeAggregatedTimeLikelihood(TimedSequence ts) {
		final TDoubleList list = computeTimeLikelihoods(ts);
		return aggregate(list, aggType);
	}

	@Override
	public boolean isAnomaly(TimedSequence s) {
		final TDoubleList eventLikelihoods = computeEventLikelihoods(s);
		final TDoubleList timeLikelihoods = computeTimeLikelihoods(s);
		return decide(eventLikelihoods, timeLikelihoods);
	}

	/**
	 * Decides whether the likelihoods indicate an anomaly
	 * 
	 * @param eventLikelihoods
	 * @param timeLikelihoods
	 * @return true for anomaly, false otherwise
	 */
	protected abstract boolean decide(TDoubleList eventLikelihoods, TDoubleList timeLikelihoods);

	@Override
	public boolean[] areAnomalies(List<TimedSequence> testSequences) {
		if (GenericSmacPipeline.isDebug()) {
			final Path testLabelFile = Paths.get("testLabels.csv");
			try {
				Files.deleteIfExists(testLabelFile);
				Files.createFile(testLabelFile);
			} catch (final IOException e1) {
				logger.error("Unexpected exception occured", e1);
			}
			try (BufferedWriter bw = Files.newBufferedWriter(testLabelFile, StandardCharsets.UTF_8)) {
				for (final TimedSequence s : testSequences) {
					bw.append(s.getLabel().toString());
					bw.append('\n');
				}
			} catch (final IOException e) {
				logger.error("Unexpected exception occured", e);
			}
		}
		final boolean[] result = new boolean[testSequences.size()];
		for (int i = 0; i < testSequences.size(); i++) {
			final TimedSequence s = testSequences.get(i);
			result[i] = isAnomaly(s);
		}
		return result;
	}

	@Override
	public void setModel(Model model) {
		if (model instanceof PDTTA) {
			logger.info("Setting model to {}", model);
			this.model = (PDTTA) model;
		} else {
			throw new UnsupportedOperationException("This AnomalyDetector can only use PDTTAs for anomaly detection");
		}
	}


	private double computeAggregatedEventLikelihood(TimedSequence s) {
		final TDoubleList list = computeEventLikelihoods(s);
		return aggregate(list, aggType);
	}

	/**
	 * 
	 * @param events
	 * @param aggType
	 * @return the list up to the last probability that exists. list may be shorter than the events list, iff there is an event which has no transition
	 */
	protected TDoubleList computeEventLikelihoods(TimedSequence s) {

		final TDoubleList list = new TDoubleArrayList();
		int currentState = 0;
		for (int i = 0; i < s.getEvents().size(); i++) {
			final Transition t = model.getTransition(currentState, s.getEvents().get(i));
			// DONE this is crap, isnt it? why not return an empty list or null iff there is no transition for the given sequence? or at least put a '0' in the
			// last slot.
			if (t == null) {
				list.add(0);
				return list;
			}
			list.add(t.getProbability());
			currentState = t.getToState();
		}
		list.add(model.getFinalStateProbability(currentState));
		return list;
	}

	public static double aggregate(TDoubleList list, ProbabilityAggregationMethod aggType) {
		double result = -1;
		if (aggType == ProbabilityAggregationMethod.MULTIPLY) {
			result = 0;
			for (int i = 0; i < list.size(); i++) {
				result += Math.log(list.get(i));
			}
		} else if (aggType == ProbabilityAggregationMethod.LUK_T) {
			result = list.get(0);
			for (int i = 1; i < list.size(); i++) {
				result = Math.max(0, list.get(i) + result - 1);
			}
		} else if (aggType == ProbabilityAggregationMethod.LUK_STRONG_DISJUNCTION) {
			result = list.get(0);
			for (int i = 1; i < list.size(); i++) {
				result = Math.min(1, list.get(i) + result);
			}
		} else if (aggType == ProbabilityAggregationMethod.NORMALIZED_MULTIPLY) {
			result = 0;
			for (int i = 0; i < list.size(); i++) {
				result += Math.log(list.get(i));
			}
			result /= list.size();
		} else if (aggType == ProbabilityAggregationMethod.NORMALIZED_MULTIPLY_UNSTABLE) {
			result = 1;
			for (int i = 0; i < list.size(); i++) {
				result *= list.get(i);
			}
			result = Math.pow(result, 1.0 / list.size());
		}
		return result;
	}
}
