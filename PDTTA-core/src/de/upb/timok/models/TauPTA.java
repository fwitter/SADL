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
package de.upb.timok.models;

import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.IntConsumer;

import jsat.distributions.Distribution;
import jsat.distributions.DistributionSearch;
import jsat.distributions.SingleValueDistribution;
import jsat.distributions.empirical.KernelDensityEstimator;
import jsat.linear.DenseVector;
import jsat.linear.Vec;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import de.upb.timok.constants.AnomalyInsertionType;
import de.upb.timok.structure.TimedSequence;
import de.upb.timok.structure.Transition;
import de.upb.timok.structure.ZeroProbTransition;

public class TauPTA extends PDTTA {
	/**
	 * 
	 */
	private static final long serialVersionUID = -7222525536004714236L;
	transient private static Logger logger = LoggerFactory.getLogger(TauPTA.class);
	TObjectIntMap<Transition> transitionCount = new TObjectIntHashMap<>();
	TIntIntMap finalStateCount = new TIntIntHashMap();

	private AnomalyInsertionType anomalyType = AnomalyInsertionType.NONE;

	int ommitedSequenceCount = 0;
	// private static final int SEQUENCE_OMMIT_THRESHOLD = 10;
	private static final double SEQUENCE_OMMIT_THRESHOLD = 0.0001;

	private static final double MAX_TYPE_FIVE_PROBABILITY = 0.2;

	public AnomalyInsertionType getAnomalyType() {
		return anomalyType;
	}

	public void setAnomalyType(AnomalyInsertionType anomalyType) {
		this.anomalyType = anomalyType;
	}

	public TauPTA(Set<Transition> transitions, TIntDoubleMap finalStateProbabilities) {
		this.transitions = transitions;
		this.finalStateProbabilities = finalStateProbabilities;
	}

	public TauPTA() {

	}

	public TauPTA(List<TimedSequence> trainingSequences) {
		super();
		final TauPTA initialPta = new TauPTA();
		initialPta.addState(START_STATE);

		for (final TimedSequence s : trainingSequences) {
			initialPta.addEventSequence(s);
		}

		// remove transitions and ending states with less than X occurences
		final double threshold = SEQUENCE_OMMIT_THRESHOLD * trainingSequences.size();
		for (final int state : initialPta.finalStateProbabilities.keys()) {
			final List<Transition> stateTransitions = initialPta.getTransitions(state, false);
			for (final Transition t : stateTransitions) {
				if (initialPta.transitionCount.get(t) < threshold) {
					initialPta.removeTransition(t);
				}
			}
			if (initialPta.finalStateCount.get(state) < threshold) {
				initialPta.finalStateCount.put(state, 0);
			}
		}

		// compute event probabilities from counts
		for (final int state : initialPta.finalStateProbabilities.keys()) {
			final List<Transition> stateTransitions = initialPta.getTransitions(state, false);
			int occurenceCount = 0;
			for (final Transition t : stateTransitions) {
				occurenceCount += initialPta.transitionCount.get(t);
			}
			occurenceCount += initialPta.finalStateCount.get(state);
			for (final Transition t : stateTransitions) {
				t.setProbability(initialPta.transitionCount.get(t) / (double) occurenceCount);
			}
			initialPta.addFinalState(state, initialPta.finalStateCount.get(state) / (double) occurenceCount);
		}
		// now the whole stuff again but only with those sequences that are in the initialPta
		// do not remove any sequences because they should occur more often than the specified threshold
		addState(START_STATE);

		for (final TimedSequence s : trainingSequences) {
			if (initialPta.isInAutomaton(s)) {
				addEventSequence(s);
			}
		}

		// compute event probabilities from counts
		for (final int state : finalStateProbabilities.keys()) {
			final List<Transition> stateTransitions = getTransitions(state, false);
			int occurenceCount = 0;
			for (final Transition t : stateTransitions) {
				occurenceCount += transitionCount.get(t);
			}
			occurenceCount += finalStateCount.get(state);
			for (final Transition t : stateTransitions) {
				t.setProbability(transitionCount.get(t) / (double) occurenceCount);
			}
			addFinalState(state, finalStateCount.get(state) / (double) occurenceCount);
		}

		// compute time probabilities
		final Map<ZeroProbTransition, TDoubleList> timeValueBuckets = new HashMap<>();
		for (final TimedSequence s : trainingSequences) {
			if (isInAutomaton(s)) {
				int currentState = START_STATE;
				for (int i = 0; i < s.length(); i++) {
					final int nextEvent = s.getEvent(i);
					final Transition t = getTransition(currentState, nextEvent);
					if (t == null) {
						// this should never happen!
						throw new IllegalStateException("Did not get a transition, but checked before that there must be transitions for this sequence " + s);
					}
					addTimeValue(timeValueBuckets, t.getFromState(), t.getToState(), t.getSymbol(), s.getTimeValue(i));
					currentState = t.getToState();
				}
			} else {
				ommitedSequenceCount++;
			}
		}
		logger.info("OmmitedSequenceCount={} out of {} sequences at a threshold of less than {} absolute occurences.", ommitedSequenceCount,
				trainingSequences.size(), SEQUENCE_OMMIT_THRESHOLD * trainingSequences.size());
		final Map<ZeroProbTransition, Distribution> distributions = fit(timeValueBuckets);
		setTransitionDistributions(distributions);
		if (distributions.size() != transitions.size()) {
			final List<Transition> missingDistributions = new ArrayList<>();
			for (final Transition t : transitions) {
				if (distributions.get(t.toZeroProbTransition()) == null) {
					missingDistributions.add(t);
				}
			}
			System.out.println(missingDistributions);
			throw new IllegalStateException("It is not possible to more/less distributions than transitions (" + distributions.size() + "/"
					+ transitions.size() + ").");
			// compute what is missing in the distribution set
		}
	}

	private void addTimeValue(Map<ZeroProbTransition, TDoubleList> result, int currentState, int followingState, int event, double timeValue) {
		final ZeroProbTransition t = new ZeroProbTransition(currentState, followingState, event);
		final TDoubleList list = result.get(t);
		if (list == null) {
			final TDoubleList tempList = new TDoubleArrayList();
			tempList.add(timeValue);
			result.put(t, tempList);
		} else {
			list.add(timeValue);
		}
	}

	private Map<ZeroProbTransition, Distribution> fit(Map<ZeroProbTransition, TDoubleList> timeValueBuckets) {
		final Map<ZeroProbTransition, Distribution> result = new HashMap<>();
		logger.debug("timevalueBuckets.size={}", timeValueBuckets.size());
		for (final ZeroProbTransition t : timeValueBuckets.keySet()) {
			result.put(t, fitDistribution(timeValueBuckets.get(t)));
		}
		return result;
	}

	@SuppressWarnings("boxing")
	private Distribution fitDistribution(TDoubleList transitionTimes) {
		final Vec v = new DenseVector(transitionTimes.toArray());
		final jsat.utils.Pair<Boolean, Double> sameValues = DistributionSearch.checkForDifferentValues(v);
		if (sameValues.getFirstItem()) {
			final Distribution d = new SingleValueDistribution(sameValues.getSecondItem());
			return d;
		} else {
			final KernelDensityEstimator kde = new KernelDensityEstimator(v);
			return kde;
		}
	}

	private void addEventSequence(TimedSequence s) {
		int currentState = START_STATE;

		for (int i = 0; i < s.length(); i++) {
			final int nextEvent = s.getEvent(i);
			Transition t = getTransition(currentState, nextEvent);
			if (t == null) {
				t = addTransition(currentState, getStateCount(), nextEvent, NO_TRANSITION_PROBABILITY);
				transitionCount.put(t, 0);
			}
			transitionCount.increment(t);
			currentState = t.getToState();
		}
		// add final state count
		finalStateCount.adjustOrPutValue(currentState, 1, 1);
	}

	// now change the pta to generate anomalies of type 1-4
	// type 1: Auf jeder Ebene des Baumes: Wähle einen zufälligen Zustand und ändere bei einer zufälligen Ausgangstransition das Symbol in ein zufälliges
	// anderes
	// Sysmbol des Alphabets für das keine andere Ausgangstranition gibt.
	// type 2: Unwahrscheinlichste Sequenzen aus PTA auswählen. nach Wahrscheinlichkeiten aller Sequenzen sortieren und die $k$ unwahrscheinlichsten Sequenzen
	// als Anomalien labeln (bzw. die Transitionen auf dem Weg der Sequenzen).
	// type 3: auf jeder Ebene des Baumes: heavily increase or decrease the outcome of one single PDF. 50%
	// type 4: Auf wahrscheinlichen Sequenzen des PTA (damit Anomalien von Typ 2 und 4 sich nicht überlappen): slightly increase or decrease (also mixed!)
	// the outcome of ALL values. 10%
	// type 5: increase or create random stop transitions? Do not increase, because it is not detectable

	// Event-Rauschen entfernen
	// now change the pta to generate anomalies of type 1-4
	// type 1: Auf jeder Ebene des Baumes: Wähle einen zufälligen Zustand und ändere bei einer zufälligen Ausgangstransition das Symbol in ein zufälliges
	// anderes
	// Sysmbol des Alphabets für das keine andere Ausgangstranition gibt.
	// type 2: Unwahrscheinlichste Sequenzen aus PTA auswählen. (Alle Sequenzen nach Wahrschienlichkeiten sortieren und die k unwahrscheinlichsten als
	// Anomalie labeln.)
	// type 3: auf jeder Ebene des Baumes: heavily increase or decrease the outcome of one single PDF. 50%
	// type 4: Auf wahrscheinlichen Sequenzen des PTA (damit Anomalien von Typ 2 und 4 sich nicht überlappen): slightly increase or decrease (also mixed!)
	// the outcome of ALL values. 10%

	public void makeAbnormal(AnomalyInsertionType newAnomalyType) {
		if (this.anomalyType != AnomalyInsertionType.NONE) {
			logger.error("A TauPTA can only have one type of anomaly. This one already has AnomalyInsertionType {}, which should be overridden with {}",
					this.anomalyType, anomalyType);
			return;
		}
		this.anomalyType = newAnomalyType;
		if (anomalyType == AnomalyInsertionType.TYPE_ONE) {
			insertPerLevelAnomaly(this::changeTransitionEvent);
		} else if (anomalyType == AnomalyInsertionType.TYPE_TWO) {
			// TODO Anomalies for type 2
			insertSequentialAnomaly(this::insertAnomaly2);
		} else if (anomalyType == AnomalyInsertionType.TYPE_THREE) {
			insertPerLevelAnomaly(this::changeTimeProbability);
		} else if (anomalyType == AnomalyInsertionType.TYPE_FOUR) {
			// TODO Anomalies for type 4
			insertSequentialAnomaly(this::insertAnomaly4);
		} else if (anomalyType == AnomalyInsertionType.TYPE_FIVE) {
			insertPerLevelAnomaly(this::addFinalStateProbability);
		} else {
			throw new IllegalArgumentException("the AnomalyInsertionType " + newAnomalyType + " is not supported.");
		}
		this.restoreConsistency();

	}

	private void insertSequentialAnomaly(Consumer consumer) {
		final List<TimedSequence> allSequences = getAllSequences();
		final TObjectDoubleMap sequenceProbabilities = computeEventProbabilities(allSequences);
		allSequences.sort((t1, t2) -> Double.compare(sequenceProbabilities.get(t1), sequenceProbabilities.get(t2)));
	}

	private TObjectDoubleMap computeEventProbabilities(List<TimedSequence> allSequences) {
		final TObjectDoubleMap result = new TObjectDoubleHashMap();
		for (final TimedSequence timedSequence : allSequences) {
			final TIntList events = timedSequence.getEvents();
			int currentState = getStartState();
			final TDoubleList probabilities = new TDoubleArrayList(events.size());
			for (int i = 0; i < events.size(); i++) {
				final int event = events.get(i);
				final Transition t = getTransition(currentState, event);
				final double probability = t.getProbability();
				probabilities.add(probability);
				currentState = t.getToState();
			}
			result.put(timedSequence, probabilities.sum());
		}
		return result;
	}

	private List<TimedSequence> getAllSequences() {
		// TODO Auto-generated method stub
		return null;
	}

	private List<TimedSequence> getAllSequences(int fromState) {
		// TODO Auto-generated method stub
		return null;
	}

	private void insertAnomaly4(Object input) {

	}

	private void insertAnomaly2(Object input) {
		// TODO Auto-generated method stub
	}




	private void insertPerLevelAnomaly(IntConsumer consumer) {
		for (int height = 0; height < getTreeHeight(); height++) {
			final TIntList states = getStates(height);
			final int chosenState = states.get(r.nextInt(states.size()));
			consumer.accept(chosenState);
		}
	}

	public static <T> T chooseRandomObject(List<T> list, Random rnd) {
		return list.get(rnd.nextInt(list.size()));
	}


	private void changeTimeProbability(int chosenState) {
		final List<Transition> possibleTransitions = getTransitions(chosenState, false);

		final Transition chosenTransition = chooseRandomObject(possibleTransitions, r);
		removeTransition(chosenTransition);

		addAbnormalTransition(chosenTransition.getFromState(), chosenTransition.getToState(), chosenTransition.getSymbol(), chosenTransition.getProbability(),
				AnomalyInsertionType.TYPE_THREE);
	}

	private void addFinalStateProbability(int chosenState) {
		// only add if there was no final state transition before
		// restore probability sum afterwards
		final List<Transition> possibleTransitions = getTransitions(chosenState, true);
		// only do so if there is no stopping transition in the possibleTransitions
		if (!possibleTransitions.stream().anyMatch(t -> t.isStopTraversingTransition())) {
			final double probability = r.nextDouble() * MAX_TYPE_FIVE_PROBABILITY;
			addFinalState(chosenState, probability);
			// now fix probs that they sum up to one
			fixProbability(chosenState);
		}
	}

	private void changeTransitionEvent(int chosenState) {
		final List<Transition> possibleTransitions = getTransitions(chosenState, false);
		final TIntList notOccuringEvents = new TIntArrayList(alphabet);
		for (final Transition t : possibleTransitions) {
			notOccuringEvents.remove(t.getSymbol());
		}
		if (notOccuringEvents.size() == 0) {
			logger.warn("Not possible to change an event in state {}", chosenState);
			return;
		}
		final Transition chosenTransition = chooseRandomObject(possibleTransitions, r);

		final int chosenEvent = notOccuringEvents.get(r.nextInt(notOccuringEvents.size()));
		removeTransition(chosenTransition);
		addAbnormalTransition(chosenTransition.getFromState(), chosenTransition.getToState(), chosenEvent, chosenTransition.getProbability(),
				AnomalyInsertionType.TYPE_ONE);
	}


	/**
	 * returns the maximum tree height
	 * 
	 * @return
	 */
	public int getTreeHeight() {
		return getTreeHeight(0, 0);
	}

	private int getTreeHeight(int currentState, int currentDepth) {
		final TIntList result = new TIntArrayList();
		final List<Transition> deeperTransitions = getTransitions(currentState, false);
		if (deeperTransitions.size() == 0) {
			return currentDepth;
		}
		for (final Transition t : deeperTransitions) {
			result.add(getTreeHeight(t.getToState(), currentDepth + 1));
		}
		return result.max();

	}

	/**
	 * returns all the states on the given tree height / tree level
	 * 
	 * @param treeHeight
	 * @return
	 */
	public TIntList getStates(int treeHeight) {
		return getStates(treeHeight, 0, 0);
	}

	private TIntList getStates(int treeHeight, int currentState, int currentDepth) {
		final TIntList result = new TIntArrayList();
		if (currentDepth < treeHeight) {
			final List<Transition> deeperTransitions = getTransitions(currentState, false);
			for (final Transition t : deeperTransitions) {
				result.addAll(getStates(treeHeight, t.getToState(), currentDepth + 1));
			}
		} else {
			result.add(currentState);
		}
		return result;
	}

}
