package sadl.models;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import sadl.constants.AnomalyInsertionType;
import sadl.input.TimedInput;
import sadl.utils.IoUtils;
import sadl.utils.MasterSeed;

public class TauPtaTestSmallV1 {
	// XXX if this test fails and V2 passes then this may be an issue with the reset method in this class. As long as V2 passes, everything should be OK
	static TimedInput trainingTimedSequences;

	@BeforeClass
	public static void setup() throws URISyntaxException, IOException {
		final Path p = Paths.get(TauPtaTestSmallV1.class.getResource("/taupta/small/rti_small.txt").toURI());
		trainingTimedSequences = TimedInput.parseAlt(p, 1);
	}

	@Before
	public void reset() {
		MasterSeed.reset();
	}



	@Test
	public void testTauPTATimedInputNormalSmall() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_normal.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}



	@Test
	public void testTauPTATimedInputAbnormal1Small() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_abnormal_1.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		pta.makeAbnormal(AnomalyInsertionType.TYPE_ONE);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}



	@Test
	public void testTauPTATimedInputAbnormal2Small() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_abnormal_2.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		pta.makeAbnormal(AnomalyInsertionType.TYPE_TWO);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}



	@Test
	public void testTauPTATimedInputAbnormal3Small() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_abnormal_3.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		pta.makeAbnormal(AnomalyInsertionType.TYPE_THREE);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}

	@Test
	public void testTauPTATimedInputAbnormal4Small() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_abnormal_4.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		pta.makeAbnormal(AnomalyInsertionType.TYPE_FOUR);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}



	@Test
	public void testTauPTATimedInputAbnormal5Small() throws IOException, URISyntaxException {
		final Path p = Paths.get(this.getClass().getResource("/taupta/small/pta_abnormal_5.xml").toURI());
		final TauPTA pta = new TauPTA(trainingTimedSequences);
		pta.makeAbnormal(AnomalyInsertionType.TYPE_FIVE);
		final TauPTA saved = (TauPTA) IoUtils.xmlDeserialize(p);
		assertEquals(pta, saved);
	}

}
