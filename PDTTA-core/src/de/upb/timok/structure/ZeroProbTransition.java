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
package de.upb.timok.structure;

import java.io.Serializable;

public class ZeroProbTransition extends Transition implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5201104088874076824L;

	public ZeroProbTransition(int fromState, int toState, int symbol) {
		super();
		this.fromState = fromState;
		this.toState = toState;
		this.symbol = symbol;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Transition other = (Transition) obj;
		if (fromState != other.fromState)
			return false;
		if (symbol != other.symbol)
			return false;
		if (toState != other.toState)
			return false;
		return true;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + fromState;
		result = prime * result + symbol;
		result = prime * result + toState;
		return result;
	}

	@Override
	public double getProbability() {
		throw new UnsupportedOperationException("This method is not supported for this class");
	}

}
