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
package de.upb.timok.constants;

import treba.trebaConstants;

public enum MergeTest {
	ALERGIA(trebaConstants.MERGE_TEST_ALERGIA), CHI_SQUARED(trebaConstants.MERGE_TEST_CHISQUARED), LR(trebaConstants.MERGE_TEST_LR), BINOMIAL(
			trebaConstants.MERGE_TEST_BINOMIAL), EXACT_M(trebaConstants.MERGE_TEST_EXACT_M), EXACT_B(trebaConstants.MERGE_TEST_EXACT_B), MDI(-1);

	private final int algorithm;

	public int getAlgorithm() {
		return algorithm;
	}

	MergeTest(int algorithm) {
		this.algorithm = algorithm;
	}
}
