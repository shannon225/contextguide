package org.searlelab.context.io;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.searlelab.context.mprophet.IsolationWindow;
import org.searlelab.msrawjava.io.encyclopedia.EncyclopeDIAFile;
import org.searlelab.context.io.TargetedBootstrapper;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;

public class TargetedBootstrapperCLI {
	
	public static void main(String[] args) throws Throwable {
		if (args.length < 3 || args.length > 7) {
			System.err.println("Usage: " + "java edu.washintgon.gs.maccoss.encyclopedia.context.TargetedBootstrapper "
					+ "<library file location> <.dia file location> <target_decoy map output location> "
					+ "\n[seed] [numberOfpeptides] [halfWindowWidthRT] [halfWindowWidthMz]");
			System.exit(1);
		}

		String libraryPath = args[0];
		String rawFilePath = args[1];
		Path mapOutputPath = Paths.get(args[2]);

		Path rawFile = Paths.get(rawFilePath);
		String baseName = rawFilePath.replaceFirst("\\.dia$", "");

		int seed = 0;
		AminoAcidConstants aaConstants = new AminoAcidConstants();
		int numberOfPeptides = 100; // number of Peptides per assay
		float halfWindowWidthRT = 0.25f;
		double halfWindowWidthMz = 1.0;

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

		for (int i = 0; i <= seed; i++) {// Randomly Select Precursors, then use them to mask the .DIA file

			ArrayList<IsolationWindow> isolationWindows = TargetedBootstrapper.selectMask(numberOfPeptides, aaConstants, i, libraryPath,
					mapOutputPath, halfWindowWidthRT);
			Path outputPath = rawFile.getParent().resolve(baseName + "_masked" + i + "_assay.dia");
			Path maskedAssayOutputPath = rawFile.getParent().resolve(baseName + "_masked" + i + "_assay.txt");

			EncyclopeDIAFile maskedFile = TargetedBootstrapper.writeMaskedFile(isolationWindows, i, rawFilePath, outputPath, halfWindowWidthMz);
			TargetedBootstrapper.writeAssayList(isolationWindows, maskedAssayOutputPath);

			System.out.println("Complete! The masked file " + maskedFile + i + " was made.\n");
		}
	}
}
