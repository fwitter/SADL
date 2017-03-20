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
package sadl.modellearner.rtiplus;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.util.Precision;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import gnu.trove.list.array.TDoubleArrayList;

public class StatisticsUtilTest {

	private static final Logger logger = LoggerFactory.getLogger(StatisticsUtilTest.class);

	private static final Random random = new Random(42);
	private static final int numTests = 1000;
	private static final int maxNumValues = 5;
	private static final double maxValue = 10.0;

	private static final double[][] valueLists = new double[numTests][];
	private static final double[][] weightLists = new double[numTests][];

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		for (int i = 0; i < numTests; i++) {
			final int length = random.nextInt(maxNumValues - 1) + 1;
			valueLists[i] = random.doubles(length, -maxValue, maxValue).map(x -> Math.round(x)).toArray();
			weightLists[i] = random.doubles(length, 0.0, maxValue).toArray();
		}
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testWeightedMedian() {

		logger.info("Starting testWeightedMedian...");

		for (int i = 0; i < numTests; i++) {
			final double wmSlow = StatisticsUtil.weightedMedianSlow(valueLists[i], weightLists[i]);
			final double wm = StatisticsUtil.weightedMedian(valueLists[i], weightLists[i]);
			Assert.assertTrue(wmSlow + " != " + wm, Precision.equals(wmSlow, wm));
		}

		logger.info("Finished testWeightedMedian.");
	}

	@Test
	public void testMedcouple() {

		logger.info("Starting testMedcouple...");

		for (int i = 0; i < numTests; i++) {
			final double[] x = Arrays.copyOf(valueLists[i], valueLists[i].length);
			Arrays.sort(x);
			System.out.println(Arrays.toString(x));
			final double median = StatisticsUtil.calculateMedian(new TDoubleArrayList(valueLists[i]), true);
			System.out.println("Median: " + median);
			final double mcSlow = StatisticsUtil.calculateMedcoupleSlow(new TDoubleArrayList(valueLists[i]), median, true);
			System.out.println("Medcouple slow: " + mcSlow);
			final double mc = StatisticsUtil.calculateMedcouple(new TDoubleArrayList(valueLists[i]), median, true);
			System.out.println("Medcouple: " + mc);
			Assert.assertTrue(mcSlow + " != " + mc, Precision.equals(mcSlow, mc));
		}

		logger.info("Finished testMedcouple.");
	}

}
