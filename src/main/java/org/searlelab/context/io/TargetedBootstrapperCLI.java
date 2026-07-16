package org.searlelab.context.io;

import java.nio.file.Path;
import java.nio.file.Paths;

public class TargetedBootstrapperCLI {
	private static final int DEFAULT_MAX_SEED = 0;
	private static final int DEFAULT_NUMBER_OF_PEPTIDES = 100;
	private static final float DEFAULT_WINDOW_WIDTH_RT = 0.50f;
	private static final double DEFAULT_HALF_WINDOW_WIDTH_MZ = 1.0;

	public static void main(String[] args) throws Throwable {
		if (args.length < 3 || args.length > 7) {
			System.err.println("Usage: " + "org.searlelab.context.io.TargetedBootstrapper "
					+ "<library file location> <.dia file location> <target_decoy map output location> "
					+ "\n[seed] [numberOfpeptides] [halfWindowWidthRT] [halfWindowWidthMz]");
			System.exit(1);
		}

		String libraryPath = args[0];
		String rawFilePath = args[1];
		Path mapOutputPath = Paths.get(args[2]);

		// Path rawFile = Paths.get(rawFilePath);
		// String baseName = rawFilePath.replaceFirst("\\.dia$", "");

		// Default parameters
		int seed = DEFAULT_MAX_SEED;
		int numberOfPeptides = DEFAULT_NUMBER_OF_PEPTIDES; // number of Peptides per assay
		float halfWindowWidthRT = DEFAULT_WINDOW_WIDTH_RT;
		double halfWindowWidthMz = DEFAULT_HALF_WINDOW_WIDTH_MZ;

		// Parameters as input
		if (args.length >= 4) {
			seed = Integer.parseInt(args[3]);
		}

		if (args.length >= 5) {
			numberOfPeptides = Integer.parseInt(args[4]);
		}

		if (args.length >= 6) {
			halfWindowWidthRT = Float.parseFloat(args[5]);
		}

		if (args.length >= 7) {
			halfWindowWidthMz = Double.parseDouble(args[6]);
		}

		TargetedBootstrapper bootstrapper = new TargetedBootstrapper();
		bootstrapper.execute(libraryPath, rawFilePath, mapOutputPath, seed, numberOfPeptides, halfWindowWidthRT,
				halfWindowWidthMz);

	}

}
