package org.searlelab.context.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import org.searlelab.context.mprophet.IsolationWindow;
import org.searlelab.msrawjava.io.encyclopedia.EncyclopeDIAFile;
import org.searlelab.msrawjava.model.FragmentScan;
import org.searlelab.msrawjava.model.PrecursorScan;
import org.searlelab.msrawjava.model.Range;
import org.searlelab.msrawjava.model.WindowData;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class TargetedBootstrapper {

	// First function - Randomly Selects Precursors from a library and compiles them
	// into a list

	public static ArrayList<IsolationWindow> selectMask(int numberOfPeptides, AminoAcidConstants aaConstants, int i,
			String libraryPath, Path mapOutputPath, float halfWindowWidthRT)
			throws IOException, SQLException, Throwable {

		// START TIMER 1
		long startTime = System.nanoTime();
		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);

		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();

		// For mapping targets and decoys later
		HashMap<String, String> targetDecoyOriginMap = new HashMap<>();

		// Set parameters before the loop
		int randomValue = 0 + i; // Add haliburton's number to get a random number

		HashSet<Integer> simulatedAssaySet = new HashSet<>();
		HashSet<String> sequencesSelectedForMasking = new HashSet<>();
		AminoAcidConstants constants = new AminoAcidConstants();
		SearchParameters params = PecanParameterParser.getDefaultParametersObject(); // need parameters to run
		// smartDecoy

		library.openFile(file);
		// Randomly select precursors loop
		try {

			// Load all entries
			ArrayList<LibraryEntry> entries = library.getAllEntries(false, aaConstants);

			while (simulatedAssaySet.size() < numberOfPeptides) {
				randomValue = RandomGenerator.randomInt(randomValue);
				int index = Math.abs(randomValue) % entries.size();
				simulatedAssaySet.add(index);
			}
			// System.out.println("Selecting " + simulatedAssaySet.size() + " precursors for
			// a fake assay.");

			for (Integer index : simulatedAssaySet) {

				// Retrieve the library entry at the random index
				LibraryEntry entry = entries.get(index);

				// Get the m/z, RT and sequence
				double targetMz = entry.getPrecursorMZ();
				float rtCenter = entry.getRetentionTimeInSec();
				String sequence = entry.getPeptideModSeq();
				byte charge = entry.getPrecursorCharge();

				// Calculate a RT ranges for the isolationWindows object
				float rtMin = (float) (rtCenter - (60 * (halfWindowWidthRT / 2)));
				float rtMax = (float) (rtCenter + (60 * (halfWindowWidthRT / 2)));

				// Add sequences to the isolationWindows object
				IsolationWindow window = new IsolationWindow(targetMz, charge, rtMin, false);
				isolationWindows.add(window);
				sequencesSelectedForMasking.add(sequence);

				// Make decoy sequences for each target peptide
				String decoy = PeptideUtils.getSmartDecoy(sequence, charge, sequencesSelectedForMasking, params);
				String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
				double decoyMz = constants.getChargedMass(correctedDecoyMass, charge);
				targetDecoyOriginMap.put(sequence, decoy);
				writeTargetDecoyMap(targetDecoyOriginMap, mapOutputPath);

				// Add decoys to Isolation Windows
				IsolationWindow decoyWindow = new IsolationWindow(decoyMz, charge, rtMin, true);
				isolationWindows.add(decoyWindow);
			}

			library.close();
			System.out.println(isolationWindows.size() + " Precursors marked for extraction.");

		} catch (Exception e) {
			System.out.println("There was an error with selecting precursors. Check file path.");
			e.printStackTrace();
		}

		// END TIMER 1
		long endTime = System.nanoTime();
		long duration = endTime - startTime;
		System.out.println("randomlySelectPrecursors(): Time taken (ms) : " + duration / 1_000_000);

		return isolationWindows;
	}

	// Second function - Uses the IsolationWindow List to mask the raw data
	@SuppressWarnings("unused")
	public static EncyclopeDIAFile writeMaskedFile(ArrayList<IsolationWindow> isolationWindows, int i, String diaFilePath,
			Path outputPath, double halfWindowWidthMz) throws Throwable {

		// START TIMER 2
		long startTime = System.nanoTime();

		File rawFile = new File(diaFilePath);
		EncyclopeDIAFile maskedFile = new EncyclopeDIAFile();
		EncyclopeDIAFile rawLibraryFile = new EncyclopeDIAFile();
		File outputFile = outputPath.toFile();

		HashSet<Integer> addedPrecursors = new HashSet<>();
		HashSet<Integer> addedFragments = new HashSet<>();

		// System.out.println("Is the .dia file open? " + rawLibraryFile.isOpen());

		try {
			rawLibraryFile.openFile(rawFile);
			maskedFile.openFile();

			// Add Ranges
			HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();
			System.out.println("Masking DIA file based on the selected precursors...");
			for (IsolationWindow window : isolationWindows) {

				double windowMz = window.getTargetMz();
				float windowStartTime = window.getRtMin();
				float windowStopTime = window.getRtMax();
				boolean sqrt = false;
				double mzStart = windowMz - halfWindowWidthMz;
				double mzStop = windowMz + halfWindowWidthMz;
				Range mzRange = new Range(mzStart, mzStop);

				ArrayList<org.searlelab.msrawjava.model.FragmentScan> fragmentScansFromWindow = rawLibraryFile.getStripes(windowMz, windowStartTime,
						windowStopTime, sqrt);
				ArrayList<FragmentScan> matchingScans = new ArrayList<>();

				// Add Fragment Scans
				for (FragmentScan scan : fragmentScansFromWindow) {
					double scanMz = scan.getPrecursorMZ();
					float scanRT = scan.getScanStartTime();
					int scanIndex = scan.getSpectrumIndex();
					if (mzRange.contains(scanMz) && !addedFragments.contains(scanIndex)) {
						matchingScans.add(scan);
						addedFragments.add(scanIndex);
					}
				}

				for (Entry<Range, WindowData> entry : rawLibraryFile.getRanges().entrySet()) {
					if (mzRange.contains(entry.getKey().getMiddle())) {
						dutyCycleMap.put(entry.getKey(), entry.getValue());
					}
				}
				maskedFile.setRanges(dutyCycleMap);
				maskedFile.addStripe(matchingScans);

				// Add Precursor Scans
				ArrayList<PrecursorScan> precursorScanFromWindow = rawLibraryFile.getPrecursors(windowStartTime,
						windowStopTime);
				ArrayList<PrecursorScan> matchingPrecursors = new ArrayList<>();

				for (PrecursorScan precursor : precursorScanFromWindow) {
					Range precursorRange = new Range(precursor.getIsolationWindowLower(),
							precursor.getIsolationWindowUpper());
					int spectrumIndex = precursor.getSpectrumIndex();
					if ((precursorRange.contains(mzRange) && !addedPrecursors.contains(spectrumIndex))) {
						matchingPrecursors.add(precursor);
						addedPrecursors.add(spectrumIndex);
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

		System.out.println("maskDIAFileBasedOnIsolationWindows(): Time taken (ms) : " + duration / 1_000_000);

		maskedFile.saveAsFile(outputFile);
		System.out.println("Target mass list for the masked file  was written to " + outputPath
				+ "\n Number of added Precursor scans: " + addedPrecursors.size()
				+ "\n Number of added Fragment scans: " + addedFragments.size());

		return maskedFile;
	}

	public static void writeAssayList(ArrayList<IsolationWindow> isolationWindows, Path outputPath) throws IOException {
		try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
			writer.write("Compound\tFormula\tAdduct\tm/z\tz\tRT Time (min)\tWindow (min)\tisDecoy");
			writer.newLine();
			for (IsolationWindow window : isolationWindows) {
				String compound = window.getCompound();
				double targetMz = window.getTargetMz();
				byte charge = window.getCharge();
				boolean isDecoy = window.isDecoy();

				float rtCenterMin = ((window.getRtMin() + window.getRtMax()) / 2.0f) / 60.0f;
				float windowMin = (window.getRtMax() - window.getRtMin()) / 60.0f;

				writer.write(compound + "\t" + "s" + "\t" + "(no adduct)" + "\t" + targetMz + "\t" + charge + "\t"
						+ rtCenterMin + "\t" + windowMin + "\t" + isDecoy);
				writer.newLine();

			}
		}
	}

	private static void writeTargetDecoyMap(HashMap<String, String> targetDecoyMap, Path outputPath)
			throws IOException {

		try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {

			writer.write("decoySequence\ttargetSequence");
			writer.newLine();

			for (Entry<String, String> entry : targetDecoyMap.entrySet()) {
				String decoySequence = entry.getKey();
				String targetSequence = entry.getValue();

				writer.write(decoySequence + "\t" + targetSequence);
				writer.newLine();
			}
		}
	}
}
