package searlelab.encyclopediana.context;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class TargetedBootstrapper {

	private static String libraryPath = "C:/Users/m334793/Documents/Library/easyspray_lit_immune_library.elib";
	private static String diaFilePath = "C:/Users/m334793/Documents/Library/cd14_combined.dia"; // Path to a new .dia file
	private static int seed = 0;
//	private static int i;
//	private static String outputPath = "C:/Users/m334793/Documents/Library/masked" + i + "_cd14_combined.dia"; // Output to new .dia file

	// @SuppressWarnings("unused")
	public static void main(String[] args) throws Throwable {
//		if (args.length!=5) {
//			Logger.errorLine("TargetedBootstrapper requires five parameters in order:");
//			Logger.logLine(" 1) Input Library (.elib file)");
//			Logger.logLine(" 2) Input raw file (.raw, .mzML, or .dia");
//			Logger.logLine(" 3) Output file path");
//			Logger.logLine(" 4) Number of peptides per assay (default = 100)");
//			Logger.logLine(" 5) Starting Seed Number (default is 3, which will use seeds 0, 1, 2, and 3, producing 4 assays total);");
		
		
//		File libraryPath = new File(args[0]);
//		File rawFile = new File(args[1]);
//		File outputPath = new File(args[2]);
//		int numberOfPeptides = Integer.parseInt(args[4]);
//		int seed = Integer.parseInt(args[5]);

		AminoAcidConstants aaConstants = new AminoAcidConstants();
		int numberOfPeptides = 100; // number of Peptides per assay


		for (int i = 0; i <= seed; i++) {// Randomly Select Precursors, then use them to mask the .DIA file

			ArrayList<IsolationWindows> isolationWindows = selectMask(numberOfPeptides, aaConstants, i);

			StripeFile maskedFile = writeMaskedFile(isolationWindows, i);

			System.out.println("Complete! The masked file " + maskedFile + i + " was made.\n");
			// System.out.println("New PRM file found at " + outputPath);
		}
	}


	// First function - Randomly Selects Precursors from a library and compiles them
	// into a list

	public static ArrayList<IsolationWindows> selectMask(int numberOfPeptides,
			AminoAcidConstants aaConstants, int i) throws IOException, SQLException, Throwable {

		// START TIMER 1
		long startTime = System.nanoTime();
		//		String libraryPath = "C:/Users/m334793/Documents/Library/easyspray_lit_immune_library.elib";

		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);
		
		ArrayList<IsolationWindows> isolationWindows = new ArrayList<>();
		int randomValue = 0 + i; // Add haliburton's number to get random number
		
		HashSet<Integer> simulatedAssaySet = new HashSet<>();
		HashSet<String> sequencesSelectedForMasking = new HashSet<>();
		AminoAcidConstants constants = new AminoAcidConstants();
		SearchParameters params = PecanParameterParser.getDefaultParametersObject(); // need parameters to run smartDecoy 

		library.openFile(file);
		//		int seed1 = 1;
		try {
			// Open library


			// Load all entries
			ArrayList<LibraryEntry> entries = library.getAllEntries(false, aaConstants);

			// Print size of library to confirm that it is open
			System.out
			.println("Total number of peptides in the chromatogram library: " + entries.size() + " peptides.");

			//		int seed = (int) 1;

			while (simulatedAssaySet.size() < numberOfPeptides) {
				randomValue = RandomGenerator.randomInt(randomValue);
				int index = Math.abs(randomValue) % entries.size();
				simulatedAssaySet.add(index);
			}
			System.out.println("Selecting " + simulatedAssaySet.size() + " precursors for a fake assay.");

			for (Integer index : simulatedAssaySet) {
				
				// Retrieve the library entry at the random index
				LibraryEntry entry = entries.get(index);
				
				// Get the m/z, RT and sequence
				double targetMz = entry.getPrecursorMZ();
				float rtCenter = entry.getRetentionTimeInSec();
				String sequence = entry.getPeptideModSeq();
				byte charge = entry.getPrecursorCharge();

				// Calculate a RT ranges for the isolationWindows object
				float rtMin = (float) (rtCenter - (60 * 2.5));
				float rtMax = (float) (rtCenter + (60 * 2.5));
				
				// Add sequences to the isolationWindows object
				IsolationWindows window = new IsolationWindows(targetMz, rtMin, rtMax, false);
				isolationWindows.add(window);
				sequencesSelectedForMasking.add(sequence);
				
				// Now let's make decoys for each target 
//				byte charge = PeptideUtils.getExpectedChargeState(sequence);
				String decoy = PeptideUtils.getSmartDecoy(sequence, charge, sequencesSelectedForMasking, params);
				String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
				double decoyMz = constants.getChargedMass(correctedDecoyMass, charge);
				
				// Add decoys to Isolation Windows
				IsolationWindows decoyWindow = new IsolationWindows(decoyMz, rtMin, rtMax, true);
				isolationWindows.add(decoyWindow);
				System.out.println("The peptide " + sequence + " had a decoy " + decoy + " made at m/z " + decoyMz + " has been selected. The decoy was generated!" 
						+ "\ndecoyMz = " + decoyMz 
						+ "\ntargetMz = " + targetMz);	
				
				System.out.println(entry.getPeptideModSeq());
				System.out.println(correctedDecoyMass);
				//

				//			System.out.println("Added the precursor at " + targetMz + " between " + rtMin / 60 + " and "
				//					+ rtMax / 60 + " minutes.");
			}
			
			// Generating decoys and adding them to the ArrayList<> isolationWindows
//			for (String sequence : sequencesSelectedForMasking) {
//
//				byte charge = PeptideUtils.getExpectedChargeState(sequence); 
//				String decoy = PeptideUtils.getSmartDecoy(sequence, charge, sequencesSelectedForMasking, params);
//				String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
//				double decoyMass = constants.getMass(correctedDecoyMass);
//				double decoyMz = decoyMass/charge;
//				
//				// Add decoys to isolation windows
//				IsolationWindows decoyWindow = new IsolationWindows(decoyMz, 0, 0, 0, true);
//				isolationWindows.add(decoyWindow);
//				
//				System.out.println("The peptide " + sequence + " at m/z " + decoyMz + " has been selected. The decoy was generated!" 
//						+ "\ndecoyMass = " + decoyMass + ", charge = +" + charge  
//						+ "\ndecoyMz = " + decoyMz);		
//			} 
 
			library.close();
			System.out.println(isolationWindows.size() + " Precursors marked for extraction.");

		} catch (Exception e) {
			System.out.println("There was an error with selecting precursors. Check file path.");
			e.printStackTrace(); // important for debugging
		}

		// END TIMER 1
		long endTime = System.nanoTime();
		long duration = endTime - startTime;
		// System.out.println("randomlySelectPrecursors(): Time taken (ns): " +
		// duration);
		System.out.println("randomlySelectPrecursors(): Time taken (ms) : " + duration / 1_000_000);

		return isolationWindows;
	}

	// Second function - Uses the IsolationWindows List to mask the raw data
	@SuppressWarnings("unused")
	public static StripeFile writeMaskedFile(ArrayList<IsolationWindows> isolationWindows, int i)
			throws Throwable {

		// START TIMER 2
		long startTime = System.nanoTime();

		File rawFile = new File(diaFilePath);
		StripeFile maskedFile = new StripeFile(false);
		StripeFile rawLibraryFile = new StripeFile(false);
		String outputPath = "C:/Users/m334793/Documents/Library/masked" + i + "_cd14_combined.dia"; // Output to new .dia file
		File outputFile = new File(outputPath);
		
		HashSet<Integer> addedPrecursors = new HashSet<>();
		HashSet<Integer> addedFragments = new HashSet<>();

		// System.out.println("Is the .dia file open? " + rawLibraryFile.isOpen());

		try {
			rawLibraryFile.openFile(rawFile);
			maskedFile.openFile();
			// Add Ranges
			HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();
			System.out.println("Masking DIA file based on the selected precursors...");
			for (IsolationWindows window : isolationWindows) {

				// System.out.println("Isolation window read from randomlySelectPrecursors().
				// TargetMz is " + window.getTargetMz() + " from " + window.getRtMin()/60 + " to
				// " + window.getRtMax()/60);	

				double windowMz = window.getTargetMz();
				float windowStartTime = window.getRtMin();
				float windowStopTime = window.getRtMax();
				boolean sqrt = false;
				double mzStart = windowMz - 0.35; // range specified for an ion trap stellar - extract out a 0.7 mz window
				double mzStop = windowMz + 0.35;
				Range mzRange = new Range(mzStart, mzStop);
				// Range rtInSecRange = new Range(mzStart*60, mzStop*60);

				ArrayList<FragmentScan> fragmentScansFromWindow = rawLibraryFile.getStripes(windowMz, windowStartTime, windowStopTime, sqrt);
				ArrayList<FragmentScan> matchingScans = new ArrayList<>();
				//				HashSet<Integer> addedFragmentSpectrumIndexes = new HashSet<>();

				// Add Fragment SCans
				for (FragmentScan scan : fragmentScansFromWindow) {
					double scanMz = scan.getPrecursorMZ();
					float scanRT = scan.getScanStartTime();
					int scanIndex = scan.getSpectrumIndex();
					if (mzRange.contains(scanMz) && !addedFragments.contains(scanIndex)) {
						matchingScans.add(scan);
						addedFragments.add(scanIndex);
					} 
				}

				for (Entry<Range, WindowData> entry : rawLibraryFile.getRanges().entrySet()) { //loadRanges instead? 
					if (mzRange.contains(entry.getKey().getMiddle())) {
						dutyCycleMap.put(entry.getKey(), entry.getValue());
					}
				}
				maskedFile.setRanges(dutyCycleMap);
				maskedFile.addStripe(matchingScans);
				// Add Precursor Scans
				ArrayList<PrecursorScan> precursorScanFromWindow = rawLibraryFile.getPrecursors(windowStartTime, windowStopTime);
				ArrayList<PrecursorScan> matchingPrecursors = new ArrayList<>();
				//				HashSet<Integer> addedPrecursorSpectrumIndexes = new HashSet<>();

				for (PrecursorScan precursor : precursorScanFromWindow) {
					Range precursorRange = new Range(precursor.getIsolationWindowLower(),
							precursor.getIsolationWindowUpper());
					//		System.out.println("Precursor range: " + precursorRange);
					int spectrumIndex = precursor.getSpectrumIndex();
					if ((precursorRange.contains(mzRange) && !addedPrecursors.contains(spectrumIndex))) {
						matchingPrecursors.add(precursor);
						addedPrecursors.add(spectrumIndex);

						//			System.out.println("Precursor count: " + matchingPrecursorScans.size());
					}
				}
				maskedFile.addPrecursor(matchingPrecursors);

			}
			maskedFile.setFileName(rawFile.getName(), null, rawFile.getAbsolutePath());
			maskedFile.addMetadata(diaFilePath, diaFilePath);
			rawLibraryFile.close();

		} catch (IOException e) {
			System.out.println("Unable to open raw file.");
			e.printStackTrace();
		}

		// END TIMER 2
		long endTime = System.nanoTime();
		long duration = endTime - startTime;

		// System.out.println("maskDIAFileBasedOnIsolationWindows(): Time taken (ns): "
		// + duration);
		System.out.println("maskDIAFileBasedOnIsolationWindows(): Time taken (ms) : " + duration / 1_000_000);

		maskedFile.saveAsFile(outputFile);
		System.out.println("Wrote the " + " was written to " + outputPath
				+ "\n Number of added Precursor scans: " + addedPrecursors.size()
				+ "\n Number of added Fragment scans: " + addedFragments.size());
		
		maskedFile.close();

		return maskedFile;
	}
}
