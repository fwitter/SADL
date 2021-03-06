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
package de.upb.timok.utils;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import weka.core.xml.XStream;

public class IoUtils {
	private static Logger logger = LoggerFactory.getLogger(IoUtils.class);

	public static void deleteFiles(String[] strings) throws IOException {
		for (final String fileName : strings) {
			final Path p = Paths.get(fileName);
			if (!Files.deleteIfExists(p)) {
				logger.warn("{} should have been explicitly deleted, but did not exist.", p);
			}
		}
	}



	public static Object deserialize(Path path) throws FileNotFoundException, IOException, ClassNotFoundException {
		try (FileInputStream fileIn = new FileInputStream(path.toFile()); ObjectInputStream in = new ObjectInputStream(fileIn)) {
			final Object o = in.readObject();
			return o;
		}
	}

	public static Object xmlDeserialize(Path path) {
		try {
			final String xml = new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
			return XStream.deSerialize(xml);
		} catch (final Exception e) {
			logger.error("Unexpected exception", e);
		}
		return null;
	}

	public static void xmlSerialize(Object o, Path path) {
		String xml;
		try {
			xml = XStream.serialize(o);
			Files.write(path, xml.getBytes());
		} catch (final Exception e) {
			logger.error("Unexpected exception", e);
		}
	}

	public static void serialize(Object o, Path path) throws IOException {
		try (FileOutputStream fileOut = new FileOutputStream(path.toFile()); ObjectOutputStream out = new ObjectOutputStream(fileOut)) {
			out.writeObject(o);
		}
	}

	public static void writeToFile(double[] testSample, Path classificationTestFile) throws IOException {
		try (BufferedWriter bw = Files.newBufferedWriter(classificationTestFile, StandardCharsets.UTF_8, StandardOpenOption.APPEND)) {
			bw.append(Arrays.toString(testSample).replace('[', ' ').replace(']', ' '));
			bw.append('\n');
		}

	}

	public static void writeToFile(List<double[]> testSamples, Path classificationTestFile) throws IOException {
		try (BufferedWriter bw = Files.newBufferedWriter(classificationTestFile, StandardCharsets.UTF_8, StandardOpenOption.APPEND)) {
			for (final double[] testSample : testSamples) {
				bw.append(Arrays.toString(testSample).replace('[', ' ').replace(']', ' '));
				bw.append('\n');
			}
		}

	}

}
