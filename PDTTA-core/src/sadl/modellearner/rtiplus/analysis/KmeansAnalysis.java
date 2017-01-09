/**
 * This file is part of SADL, a library for learning all sorts of (timed) automata and performing sequence-based anomaly detection.
 * Copyright (C) 2013-2017  the original author or authors.
 *
 * SADL is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * SADL is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with SADL.  If not, see <http://www.gnu.org/licenses/>.
 */
package sadl.modellearner.rtiplus.analysis;

import java.util.List;

import jsat.DataSet;
import jsat.SimpleDataSet;
import jsat.classifiers.DataPoint;
import jsat.clustering.SeedSelectionMethods.SeedSelection;
import jsat.clustering.kmeans.HamerlyKMeans;
import jsat.clustering.kmeans.KMeans;
import jsat.linear.distancemetrics.DistanceMetric;
import sadl.utils.MasterSeed;

public class KmeansAnalysis extends JsatClusteringAnalysis {

	private final SeedSelection seedSel;
	private final int numClusters;

	public KmeansAnalysis(DistanceMetric dstMetric, SeedSelection seedSelection, int numClusters, double clusterExpansionRate) {
		this(dstMetric, seedSelection, numClusters, clusterExpansionRate, null, -1);
	}

	public KmeansAnalysis(DistanceMetric dstMetric, SeedSelection seedSelection, int numClusters, double clusterExpansionRate,
			DistributionAnalysis fewElementsAnalysis, int fewElementsLimit) {
		super(dstMetric, clusterExpansionRate, fewElementsAnalysis, fewElementsLimit);

		this.seedSel = seedSelection;
		this.numClusters = numClusters;
	}

	@Override
	List<List<DataPoint>> computeJsatClusters(DataSet<SimpleDataSet> ds) {

		final KMeans k = new HamerlyKMeans(dm, seedSel, MasterSeed.nextRandom());
		return k.cluster(ds, numClusters);
	}

}
