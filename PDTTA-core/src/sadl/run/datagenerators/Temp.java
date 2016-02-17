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
package sadl.run.datagenerators;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Consumer;
import java.util.regex.Matcher;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import sadl.constants.Algoname;
import sadl.constants.AnomalyInsertionType;
import sadl.constants.EventsCreationStrategy;
import sadl.constants.KDEFormelVariant;
import sadl.input.TimedInput;
import sadl.input.TimedWord;
import sadl.modellearner.ButlaPdtaLearner;
import sadl.modellearner.TauPtaLearner;
import sadl.models.TauPTA;
import sadl.models.pta.Event;
import sadl.utils.CollectionUtils;
import sadl.utils.IoUtils;
import sadl.utils.MasterSeed;

public class Temp {

	private static Path DATA_TYPE = Paths.get("real");
	public static final String TRAIN_TEST_SEP = "?????????????????????????";
	private static double ANOMALY_PERCENTAGE = 0.1;
	private static int TRAIN_SIZE = 9000;
	private static int TEST_SIZE = 4000;
	private static final int SAMPLE_FILES = 11;
	Random r;

	private static Logger logger = LoggerFactory.getLogger(Temp.class);

	public static void main(String[] args) throws IOException {
		final Temp t = new Temp();
		t.temp(args[0]);
	}

	void temp(String dataString) throws IOException {
		final Path outputDir = Paths.get("output2");
		IoUtils.cleanDir(outputDir);
		if (r == null) {
			r = MasterSeed.nextRandom();
		}
		final Path confFile = Paths.get("conf-template.txt");
		final DecimalFormat df = new DecimalFormat("00");
		final EventsCreationStrategy[] splitEvents = { EventsCreationStrategy.DontSplitEvents, EventsCreationStrategy.SplitEvents };
		TimedInput trainingTimedSequences;
		TimedInput splitTrainingTimedSequences = null;
		logger.info("Parsing input file {}...", dataString);
		final TimedInput unsplitTrainingTimedSequences = TimedInput.tryParse(Paths.get(dataString));
		if (unsplitTrainingTimedSequences.size() > 13000) {
			DATA_TYPE = Paths.get("real");
		} else if (unsplitTrainingTimedSequences.size() >= 1000) {
			final String confSubFolder = "synthetic-pDtta";
			DATA_TYPE = Paths.get(confSubFolder);
			TRAIN_SIZE = 700;
			TEST_SIZE = 300;
			ANOMALY_PERCENTAGE = .5;
		} else {
			throw new IllegalStateException("Input must at least have size 1000, but has " + unsplitTrainingTimedSequences.size());
		}
		logger.info("Parsed input file.");
		for (final AnomalyInsertionType type : AnomalyInsertionType.values()) {
			if (type == AnomalyInsertionType.NONE) {
				continue;
			}

			final Path confDir = outputDir.resolve("confs").resolve(DATA_TYPE);
			for (final EventsCreationStrategy split : splitEvents) {
				if (split == EventsCreationStrategy.SplitEvents) {
					if (splitTrainingTimedSequences == null) {
						final ButlaPdtaLearner butla = new ButlaPdtaLearner(10000, EventsCreationStrategy.SplitEvents, KDEFormelVariant.OriginalKDE);
						logger.info("Splitting input into subevents...");
						final Pair<TimedInput, Map<String, Event>> p = butla.splitEventsInTimedSequences(unsplitTrainingTimedSequences);
						splitTrainingTimedSequences = p.getKey();
						logger.info("Split input into subevents.");
					}
					trainingTimedSequences = splitTrainingTimedSequences;
				} else {
					trainingTimedSequences = unsplitTrainingTimedSequences;
				}
				final String genFolder = "tpta-" + (split == EventsCreationStrategy.SplitEvents ? "prep" : "noPrep");
				final int typeIndex = type.getTypeIndex();
				final String typeFolderString = typeIndex == AnomalyInsertionType.ALL.getTypeIndex() ? "mixed" : "type" + typeIndex;
				final Path typeFolder = DATA_TYPE.resolve(genFolder).resolve(Paths.get(typeFolderString));
				final TimedInput input = trainingTimedSequences;
				final Consumer<Path> anomalyGenerator = (Path dataOutputFile) -> {
					try {
						generateModelAnomaly(input, dataOutputFile, type);
					} catch (final Exception e) {
						e.printStackTrace();
						logger.error("Unexpected exception", e);
					}
				};
				createFiles(outputDir, confFile, DATA_TYPE, df, confDir, genFolder, typeFolderString, typeFolder, anomalyGenerator);
			}

			// randomly
			final String genFolder = "random";
			final int typeIndex = type.getTypeIndex();
			final String typeFolderString = typeIndex == AnomalyInsertionType.ALL.getTypeIndex() ? "mixed" : "type" + typeIndex;
			final Path typeFolder = DATA_TYPE.resolve(genFolder).resolve(Paths.get(typeFolderString));
			final Consumer<Path> anomalyGenerator = (Path dataOutputFile) -> {
				try {
					generateRandomAnomaly(unsplitTrainingTimedSequences, dataOutputFile, type);
				} catch (final Exception e) {
					e.printStackTrace();
					logger.error("Unexpected exception", e);
				}
			};
			createFiles(outputDir, confFile, DATA_TYPE, df, confDir, genFolder, typeFolderString, typeFolder, anomalyGenerator);
		}
	}



	private void generateModelAnomaly(TimedInput input, Path dataOutputFile, AnomalyInsertionType type) throws IOException {
		if (type != AnomalyInsertionType.ALL) {
			generateSingleModelAnomaly(input, dataOutputFile, type);
		} else {
			generateMixedModelAnomaly(input, dataOutputFile);
		}
	}

	public void createFiles(final Path outputDir, final Path confFile, final Path dataType, final DecimalFormat df, final Path confDir,
			final String genFolder, final String typeFolderString, final Path typeFolder, Consumer<Path> anomalyGenerator) throws IOException {
		for (int k = 0; k < SAMPLE_FILES; k++) {
			final String destFolder = k == 0 ? "train" : "test";
			final Path dataFolder = outputDir.resolve("smac-data").resolve(dataType).resolve(genFolder).resolve(typeFolderString)
					.resolve(destFolder);
			final Path dataOutputFile = dataFolder.resolve(Paths.get(genFolder + "-" + df.format(k) + "_smac_" + typeFolderString + ".txt"));
			if (Files.notExists(dataFolder)) {
				Files.createDirectories(dataFolder);
			}
			anomalyGenerator.accept(dataOutputFile);
			logger.info("Wrote file #{} ({})", Integer.toString(k), dataOutputFile);
		}
		// final String fileSuffix = "-" + dataType + "-" + genFolder + "-" + typeFolderString + ".txt";
		// createConfDir(confFile, confDir, typeFolder, fileSuffix);
	}

	public void createConfDir(final Path confFile, final Path confDir, final Path typeFolder, final String fileSuffix) throws IOException {
		for (final Algoname algo : Algoname.values()) {
			if (Files.notExists(confDir.resolve(algo.name().toLowerCase()))) {
				Files.createDirectories(confDir.resolve(algo.name().toLowerCase()));
				logger.info("Created directory {}", confDir.resolve(algo.name().toLowerCase()));
			}
			final String fileName = algo.name().toLowerCase() + fileSuffix;
			final Path destFile = confDir.resolve(algo.name().toLowerCase()).resolve(fileName);
			Files.copy(confFile, destFile);
			final List<String> lines = Files.readAllLines(destFile);
			for (int i = 0; i < lines.size(); i++) {
				String line = lines.get(i);
				line = line.replaceAll("\\$algoname", Matcher.quoteReplacement(algo.name().toLowerCase()));
				line = line.replaceAll("\\$trainFolder", Matcher.quoteReplacement(typeFolder.resolve("train").toString()));
				line = line.replaceAll("\\$testFolder", Matcher.quoteReplacement(typeFolder.resolve("test").toString()));
				lines.set(i, line);
			}
			Files.write(destFile, lines, StandardOpenOption.WRITE);
			logger.info("Wrote file {}", destFile);
		}
	}

	Map<Path, Pair<TauPTA, TauPTA>> singleModelPtas = new HashMap<>();
	private void generateSingleModelAnomaly(TimedInput trainingTimedSequences, Path dataOutputFile, AnomalyInsertionType type) throws IOException {
		final List<TimedWord> trainSequences = new ArrayList<>();
		final List<TimedWord> testSequences = new ArrayList<>();
		final Path key = dataOutputFile.getParent().getParent();
		final Pair<TauPTA, TauPTA> result = singleModelPtas.get(key);
		TauPTA pta;
		TauPTA anomaly;
		if (result == null) {
			logger.info("Learning TPTA for single anomaly of type {}...", type);
			final TauPtaLearner learner = new TauPtaLearner();
			pta = learner.train(trainingTimedSequences);
			anomaly = SerializationUtils.clone(pta);
			logger.info("inserting Anomaly Type {} into tpta", type);
			anomaly.makeAbnormal(type);
			if (type == AnomalyInsertionType.TYPE_TWO) {
				anomaly.removeAbnormalSequences(pta);
			}
			for (int i = 0; i < TRAIN_SIZE; i++) {
				trainSequences.add(pta.sampleSequence());
			}
			singleModelPtas.put(key, Pair.of(pta, anomaly));
		} else {
			logger.info("Using cached TPTA.");
			pta = result.getKey();
			anomaly = result.getValue();
		}

		// PTAs of Type 2 and 4 always produce abnormal sequences
		// it is possible to sample abnormal and normal sequences with abnormal ptas of the other types (1,3,5).
		// but I don't know how the distribution is, so to be fair, i sample all anomalies the same
		for (int i = 0; i < TEST_SIZE; i++) {
			if (r.nextDouble() < ANOMALY_PERCENTAGE) {
				boolean wasAnormal = false;
				TimedWord seq = null;
				while (!wasAnormal) {
					seq = anomaly.sampleSequence();
					wasAnormal = seq.isAnomaly();
				}
				testSequences.add(seq);
			} else {
				testSequences.add(pta.sampleSequence());
			}
		}
		IoUtils.writeTrainTestFile(dataOutputFile, new TimedInput(trainSequences), new TimedInput(testSequences));
	}

	int[] intIndex;
	public void generateRandomAnomaly(TimedInput trainingTimedSequences, Path dataOutputFile, AnomalyInsertionType type) throws IOException {
		logger.info("generating random anomalies for type {}...", type);
		final ArrayList<TimedWord> trainSequences = new ArrayList<>();
		final ArrayList<TimedWord> testSequences = new ArrayList<>();
		// final Path p = Paths.get("pta_normal.dot");
		// pta.toGraphvizFile(outputDir.resolve(p), false);
		// final Process ps = Runtime.getRuntime().exec("dot -Tpdf -O " + outputDir.resolve(p));
		// System.out.println(outputDir.resolve(p));
		// ps.waitFor();
		if (intIndex == null || intIndex.length != trainingTimedSequences.size()) {
			intIndex = new int[trainingTimedSequences.size()];
			for (int i = 0; i < trainingTimedSequences.size(); i++) {
				intIndex[i] = i;
			}
		}
		final TIntList shuffledIndex = new TIntArrayList(Arrays.copyOf(intIndex, intIndex.length));
		shuffledIndex.shuffle(r);
		for (int i = 0; i < TRAIN_SIZE; i++) {
			trainSequences.add(trainingTimedSequences.get(shuffledIndex.get(i)));
		}
		for (int i = TRAIN_SIZE; i < TRAIN_SIZE + TEST_SIZE; i++) {
			testSequences.add(trainingTimedSequences.get(shuffledIndex.get(i)));
		}
		if (type != AnomalyInsertionType.NONE) {
			logger.info("inserting random  Anomaly Type {}", type);
			final List<TimedWord> trainSequenceClone = SerializationUtils.clone(trainSequences);
			final List<TimedWord> testSequenceClone = SerializationUtils.clone(testSequences);
			final TimedInput trainSet = new TimedInput(trainSequenceClone);
			TimedInput testSet = new TimedInput(testSequenceClone);
			testSet = testSet.insertRandomAnomalies(type, ANOMALY_PERCENTAGE);
			IoUtils.writeTrainTestFile(dataOutputFile, trainSet, testSet);

		}
	}

	private List<TauPTA> createAbnormalPtas(TauPTA normalPta) {
		final List<TauPTA> abnormalPtas = new ArrayList<>();
		for (final AnomalyInsertionType type : AnomalyInsertionType.values()) {
			if (type != AnomalyInsertionType.NONE && type != AnomalyInsertionType.ALL) {
				if (type == AnomalyInsertionType.TYPE_TWO) {
					normalPta = SerializationUtils.clone(normalPta);
					normalPta.setRandom(MasterSeed.nextRandom());
				}
				final TauPTA anomaly = SerializationUtils.clone(normalPta);
				logger.info("inserting Anomaly Type {}", type);
				anomaly.makeAbnormal(type);
				abnormalPtas.add(type.getTypeIndex() - 1, anomaly);
				if (type == AnomalyInsertionType.TYPE_TWO) {
					anomaly.removeAbnormalSequences(normalPta);
				}
			}
		}
		return abnormalPtas;
	}

	Map<Path, Pair<TauPTA, List<TauPTA>>> mixedModelPtas = new HashMap<>();
	private void generateMixedModelAnomaly(TimedInput input, Path dataOutputFile) throws IOException {
		final Path key = dataOutputFile.getParent().getParent();
		final Pair<TauPTA, List<TauPTA>> result = mixedModelPtas.get(key);
		TauPTA normalPta;
		final List<TauPTA> abnormalPtas;
		if (result == null) {
			logger.info("Learning TPTA for mixed anomalies...");
			final TauPtaLearner learner = new TauPtaLearner();
			normalPta = learner.train(input);
			abnormalPtas = createAbnormalPtas(normalPta);
			mixedModelPtas.put(key, Pair.of(normalPta, abnormalPtas));
			logger.info("Finished TauPTA creation.");
		} else {
			logger.info("Using cached TPTA.");
			normalPta = result.getKey();
			abnormalPtas = result.getValue();
		}
		final List<TimedWord> trainSequences = new ArrayList<>();
		final List<TimedWord> testSequences = new ArrayList<>();
		// final Path p = Paths.get("pta_normal.dot");
		// pta.toGraphvizFile(outputDir.resolve(p), false);
		// final Process ps = Runtime.getRuntime().exec("dot -Tpdf -O " + outputDir.resolve(p));
		// System.out.println(outputDir.resolve(p));
		// ps.waitFor();
		for (int i = 0; i < TRAIN_SIZE; i++) {
			trainSequences.add(normalPta.sampleSequence());
		}
		for (int i = 0; i < TEST_SIZE; i++) {
			if (r.nextDouble() < ANOMALY_PERCENTAGE) {
				boolean wasAnormal = false;
				TimedWord seq = null;
				final TauPTA chosen = CollectionUtils.chooseRandomObject(abnormalPtas, r);
				while (!wasAnormal) {
					seq = chosen.sampleSequence();
					wasAnormal = seq.isAnomaly();
				}
				testSequences.add(seq);
			} else {
				testSequences.add(normalPta.sampleSequence());
			}
		}
		IoUtils.writeTrainTestFile(dataOutputFile, new TimedInput(trainSequences), new TimedInput(testSequences));
	}
}
