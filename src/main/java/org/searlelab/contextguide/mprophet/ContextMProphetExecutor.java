package org.searlelab.contextguide.mprophet;

import java.io.File;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

public class ContextMProphetExecutor {

	public static void main(String[] args) {

		// Map files 
		String libraryPath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_1min/IL2_and_IL15_Combo.elib";
		String fastaPath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_1min/mus_musculus_reviewed_uniprot.fasta";

		// Where the feature files are located: 
		String diaFolderPath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_1min/";

		// Get a list of .dia files 
		File diaFolder = new File(diaFolderPath);
		File[] diaFiles = diaFolder.listFiles();

		System.out.println("DIA Folder was identified as: " + diaFolder.getAbsolutePath());
		System.out.println("DIA files detected! The following files will be processed with MProphet:" + diaFolder.listFiles());

		if (diaFiles != null) {
			for (File diaFile : diaFiles) {
				if (diaFile.isFile() && diaFile.getName().endsWith(".dia")) {
					System.out.println(diaFile.getName());
				}
			}
		}

		// Identify the current file so we can loop through all files

		if (diaFiles !=null) {
			for (File diaFile : diaFiles) {

				// Ignore files that do not end in .dia 
				if (!diaFile.isFile() || !diaFile.getName().endsWith(".dia")) {
					continue;
				}

				String currentDiaFilePath = diaFile.getAbsolutePath();
				String diaName = diaFile.getName(); 
				String baseName = diaName.substring(0, diaName.lastIndexOf(".dia"));

				File massListFile = new File(diaFolder, baseName + ".txt");
				String massListPath = massListFile.getAbsolutePath();

				System.out.println("Processesing " + diaFile.getName());

				if (!massListFile.exists()) {
					System.out.println("Skipping " + diaFile.getName() + " because mass list was not found: " + massListPath);
					continue;
				}
				executeContextMProphet(libraryPath, fastaPath, currentDiaFilePath, massListPath, diaFolder);
			}	
		}
	}

	public static void executeContextMProphet(String libraryPath, String fastaPath, String diaFilePath, String massListPath, File diaFolder) {
		File fasta = new File(fastaPath);
		File diaFile = new File(diaFilePath);
		File library = new File(libraryPath);

		String baseName = diaFilePath.replaceFirst("\\.dia$", "");

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();


		// Score features in the .dia file against the library, split the results
		try {
			ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
			String featureFileName = baseName.replaceAll("\\.txt$", "");

			File backgroundFeatureFile = new File(featureFileName + "_background.features.txt");
			File referenceFeatureFile = new File(featureFileName + "_reference.features.txt");

			MProphetExecutionData backgroundData = makeMProphetExecutionData(backgroundFeatureFile, fasta, params, ".pep");
			MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");

			float peptideFDRThreshold = 0.01f;
			int seed = 1;
			int round = 1;

			MProphetResult backgroundMProphetResult = MProphetReiter.executeMProphetTSV(backgroundData, peptideFDRThreshold, seed, params.getAAConstants(), round);
			LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();

			// 	Use the background LDA model on the reference feature file without retraining
			MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());

			System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + referenceFeatureFile.getAbsolutePath());
			//			System.out.println("Background passing peptides: " + backgroundMProphetResult.getPassingPeptides().size());
			System.out.println("Finished scoring peptides with background-trained lda model. "
					+ "\nReference passing peptides: " + referenceMProphetResult.getPassingPeptides().size());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}



	public static void executeMProphet(String libraryPath, String fastaPath, String diaFilePath, String massListPath, File diaFolder) {
		File fasta = new File(fastaPath);
		File diaFile = new File(diaFilePath);
		File library = new File(libraryPath);

		String baseName = diaFilePath.replaceFirst("\\.dia$", "");

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();


		// Score features in the .dia file against the library, split the results
		try {
			ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
			String featureFileName = baseName.replaceAll("\\.txt$", "");

			File featureFile = new File(featureFileName + "_mprophet.features.txt");
//			File referenceFeatureFile = new File(featureFileName + "_reference.features.txt");

			MProphetExecutionData featureData = makeMProphetExecutionData(featureFile, fasta, params, ".pep");
//			MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");

			float peptideFDRThreshold = 0.01f;
			int seed = 1;
			int round = 1;

			MProphetResult mprophetResult = MProphetReiter.executeMProphetTSV(featureData, peptideFDRThreshold, seed, params.getAAConstants(), round);
	//		LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();

			// 	Use the background LDA model on the reference feature file without retraining
	//		MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());

			System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + featureFile.getAbsolutePath());
			//			System.out.println("Background passing peptides: " + backgroundMProphetResult.getPassingPeptides().size());
			System.out.println("Finished scoring peptides with background-trained lda model. "
					+ "\nReference passing peptides: " + mprophetResult.getPassingPeptides().size());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}


	private static  MProphetExecutionData makeMProphetExecutionData(File inputFeatureFile, File fasta, SearchParameters params, String outputSuffix) {

		File peptideOutputFile = new File(inputFeatureFile.getAbsolutePath().replaceAll("\\.txt$", "") + outputSuffix + ".output.txt");
		File peptideDecoyFile = new File(inputFeatureFile.getAbsolutePath().replaceAll("\\.txt$", "") + outputSuffix + ".decoy.txt");

		return new MProphetExecutionData(
				inputFeatureFile,
				fasta,
				peptideOutputFile,
				peptideDecoyFile,
				params);
	}
}


