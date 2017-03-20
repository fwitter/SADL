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
public class NewOutlierAnalysis extends OutlierDistanceAnalysis {

	private final boolean onlyFarOuts;
	private final boolean presumeOnlyFewOutliers;

	public NewOutlierAnalysis(double strength, boolean onlyFarOuts, boolean presumeOnlyFewOutliers) {
		this(strength, onlyFarOuts, presumeOnlyFewOutliers, null, -1);
	}

	public NewOutlierAnalysis(double strength, boolean onlyFarOuts, boolean presumeOnlyFewOutliers, DistributionAnalysis fewElementsAnalysis,
			int fewElementsLimit) {
		super(strength, fewElementsAnalysis, fewElementsLimit);
		this.onlyFarOuts = onlyFarOuts;
		this.presumeOnlyFewOutliers = presumeOnlyFewOutliers;
	}

	@Override
	int getOutlierDistance(TDoubleList distValues) {

		// http://images.alfresco.advanstar.com/alfresco_images/pharma/2014/08/22/56d1d667-8be0-489b-8d56-dfd994f0b6f1/article-4509.pdf
		// http://academic.uprm.edu/eacuna/paperout.pdf
		// http://www.umiacs.umd.edu/labs/cvl/pirl/vikas/Software/optimal_bw/optimal_bw_code.htm
		// http://homepages.inf.ed.ac.uk/imurray2/pub/09gpds/gpds_nips.pdf
		// https://en.wikipedia.org/wiki/Dirichlet_process
		// http://jmlr.org/papers/volume13/kim12b/kim12b.pdf
		// http://www.ssc.wisc.edu/~bhansen/papers/ncde.pdf
		// http://www.bundesbank.de/Redaktion/EN/Downloads/Publications/Discussion_Paper_1/2013/2013_02_01_dkp_02.pdf?__blob=publicationFile
		// http://www.ise.bgu.ac.il/faculty/mlast/papers/outliers2.pdf
		// http://www.ifau.se/globalassets/pdf/se/2001/wp01-07.pdf

		return 0;
	}

}
