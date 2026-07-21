package org.searlelab.context.mprophet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

import org.searlelab.msrawjava.io.StripeFileInterface;
import org.searlelab.context.encyclopedia.EncyclopediaTwo;
import org.searlelab.context.encyclopedia.SearchToBLIB;
import org.searlelab.context.io.IsolationWindowReader;

public class ContextFeatureScorer {

	public static void parseCLIForScoringFeatures(String[] args) throws IOException, SQLException, InterruptedException, DataFormatException {

		if (args.length != 4) {
			System.err.println("Usage: ");
			System.err.println("java edu.washington.gs.maccoss.encyclopedia.context.ContextFeatureScorer " + 
					"<rawFilePath> <libraryFilePath> <fastaPath> <massListPath");
		}

		String rawFilePath = args[0];
		String libraryFilePath = args[1];
		String fastaPath = args[2];
		String massListPath = args[3];

		String baseName = rawFilePath.replaceFirst("\\.dia$", "");
		final File fasta = new File(fastaPath);
		File rawFile = new File(rawFilePath);
		File library = new File(libraryFilePath);

		try {
			ArrayList<ScoredFeature> partitionedFeatures = scoreFeatures(library, rawFile, fasta, baseName, massListPath);
			System.out.println("Features are scored and paritioned at " + partitionedFeatures.getLast());
		} catch (Exception e) {
			System.out.println("Something did not work... see the error tace");
			e.printStackTrace();
		} finally {
		System.out.println("Program concluded.");
		}
		}

	static IsolationWindow findMatchingMassListWindow(ScoredFeature feature, ArrayList<IsolationWindow> targetWindows) {
		String sequence = cleanPeptideSequence(feature.getSequence());

		for (IsolationWindow window : targetWindows) {
			String compound = cleanPeptideSequence(window.getCompound());

			if (compound.equals(sequence)) {
				return window;
			}
		}

		return null;
	}


	private static String cleanPeptideSequence(String sequence) {
		if (sequence == null) return "";

		return sequence.trim().replaceFirst("^[A-Za-z-]?\\.", "").replaceFirst("\\.[A-Za-z-]?$", "");
	}

	private static void writeScoredFeatures(File outputFile, ArrayList<ScoredFeature> features, String header)
			throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
			writer.write(header);
			writer.newLine();

			for (ScoredFeature feature : features) {
				writer.write(feature.getOriginalLine());
				writer.newLine();
			}
		}
	}

	// Changed isFeatureOnMassList to only check for peptide sequence equivalence, not for mass, charge and RT equivalence. I don't think is needed, but keeping for now. 
	private static byte parseCharge(String[] columns) {
		if (Integer.parseInt(columns[23])==1) return 1;
		if (Integer.parseInt(columns[24])==1) return 2;
		if (Integer.parseInt(columns[25])==1) return 3;
		if (Integer.parseInt(columns[26])==1) return 4;

		throw new IllegalArgumentException("Charge was unable to be detected from the processed features. Check input file.");
	}

	public static ArrayList<ScoredFeature> scoreFeatures(File library, File rawFile, File fasta, String baseName,
			String massListPath) throws IOException, SQLException, DataFormatException, InterruptedException {

		// Run an Encyclopedia job
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();
		LibraryScoringFactory scoringForLibrary = EncyclopediaScoringFactory.getDefaultScoringFactory(params);
		LibraryInterface interfaceForLibrary = BlibToLibraryConverter.getFile(library, fasta, params);

		EncyclopediaTwoJobData job = new EncyclopediaTwoJobData(rawFile, fasta, interfaceForLibrary,
				interfaceForLibrary, rawFile, scoringForLibrary);

		ProgressIndicator progress = new EmptyProgressIndicator(true);
		org.searlelab.msrawjava.io.StripeFileInterface interfaceForStripeFile = job.getDiaFileReader();

		// Run Encyclopedia job to get the feature file
		File featuresToSplit = job.getPercolatorFiles().getInputTSV();

		if (featuresToSplit.exists() && featuresToSplit.canRead()) {
			System.out.println("Feature file already exists, skipping feature calculation!");
			System.out.println(featuresToSplit.getAbsolutePath());
		} else {
			System.out.println("Calculating features...");
			EncyclopediaTwo.generateFeatureFile(progress, interfaceForLibrary, job,
					interfaceForStripeFile, java.util.Optional.empty());
		}

		ArrayList<ScoredFeature> uniqueFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> uniqueFeaturesList = uniqueFeatures;
		HashMap<String, ScoredFeature> bestFeatureByPeptide = new HashMap<>();
		String header;

		// Read all rows the feature file
		try (BufferedReader br = new BufferedReader(new FileReader(featuresToSplit))) {
			header = br.readLine();
			if (header == null) {
				throw new IOException("Feature file is empty, so no header could be read.");
			}

			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split("\t", -1);

				double mz = Double.parseDouble(columns[27]);
				byte featureCharge = parseCharge(columns);
				boolean isDecoy = Integer.parseInt(columns[1]) == -1;
				float primary = Float.parseFloat(columns[3]);
				float retentionTime = Float.parseFloat(columns[29]);
				String sequence = columns[30];
				String protein = columns[31];

				ScoredFeature feature = new ScoredFeature(mz, featureCharge, isDecoy, primary, retentionTime, sequence, protein, line);

				uniqueFeaturesList.add(feature);

				ScoredFeature currentBest = bestFeatureByPeptide.get(sequence);

				// Take the peptide with a higher primary score and place it on a new list

				if (currentBest == null || feature.getPrimary() > currentBest.getPrimary()) {
					bestFeatureByPeptide.put(sequence, feature);
				}

			}
		}
		ArrayList<ScoredFeature> bestFeatures = new ArrayList<>(bestFeatureByPeptide.values());
		bestFeatures.sort(Comparator.comparing(ScoredFeature::getPrimary).reversed());

		System.out.println("Selecting the betst feature per peptide..." + bestFeatures.size() + " peptides remain.");


		// Output Paths
		String referenceOutputPath = baseName + "_reference.features.txt";
		String backgroundOutputPath = baseName + "_background.features.txt";

		// Output Files
		File referenceOutput = new File(referenceOutputPath);
		File backgroundOutput = new File(backgroundOutputPath);

		// Target mass list
		ArrayList<IsolationWindow> targetWindows = IsolationWindowReader.parseMassList(massListPath);
		System.out.println(targetWindows.size() + " windows cataloged from the mass list");

		ArrayList<ScoredFeature> referenceFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> backgroundFeatures = new ArrayList<>();

		ArrayList<ScoredFeature> partitionedFeatures = new ArrayList<>(); // so that the return is all of the features

		for (ScoredFeature feature : bestFeatures) {
			IsolationWindow matchingWindow = findMatchingMassListWindow(feature, targetWindows);

			boolean isOnMassList = matchingWindow != null;
	//		boolean isMassListDecoy = isOnMassList && matchingWindow.isDecoy();
			boolean isBackground = !isOnMassList;

			ScoredFeature annotatedFeature = new ScoredFeature(feature.getMz(), feature.isDecoy(), feature.getPrimary(), feature.getRetentionTime(), cleanPeptideSequence(feature.getSequence()), feature.getProtein(), feature.getOriginalLine(), isBackground);
			partitionedFeatures.add(annotatedFeature);

			//			System.out.println("Features are being read, sequence for this feature is " + cleanPeptideSequence(feature.getSequence()));

			if (isOnMassList) {
				referenceFeatures.add(feature); 
			} else {
				backgroundFeatures.add(feature);
			}
		}
		writeScoredFeatures(referenceOutput, referenceFeatures, header);
		writeScoredFeatures(backgroundOutput, backgroundFeatures, header);

		System.out.println("Reference target features: " + referenceFeatures.size());
		//		System.out.println("Background target features: " + backgroundFeatures.size());
		//		System.out.println("Reference decoy features: " + referenceDecoyFeatures.size());
		//		System.out.println("Background decoy features: " + backgroundDecoyFeatures.size());
		return partitionedFeatures;
	}


}
