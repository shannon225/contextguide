package org.searlelab.context.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.context.mprophet.IsolationWindow;

public class IsolationWindowReader {

	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);
		boolean hasPrintedDebugInfo = false;
		boolean hasPrintedAddingPrecursorInfo = false; 

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
			String delim = getDelimiter(massListFile);
	
			@SuppressWarnings("unused")
			String header = br.readLine();
			
			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split(delim, -1); // identify the format of the mass list
				
				// Loop so that the message is only printed once
				if (!hasPrintedDebugInfo) {
					System.out.println("Reading in the window:  " + line); // Console will print what the data looks like as its read in

					hasPrintedDebugInfo = true;
				}
				
				String peptide = columns[0];  // read columns in as primitives
				double targetMz = Double.parseDouble(columns[3]);
				byte charge = Byte.parseByte(columns[4]);
				float rtCenter = Float.parseFloat(columns[5]);
				float rtWindow = Float.parseFloat(columns[6]);

				float rtMin = (rtCenter - (rtWindow / 2))*60;
				float rtMax = (rtCenter + (rtWindow / 2))*60;

				boolean isDecoy = Boolean.parseBoolean(columns[7]);

				// Assemble each window
				IsolationWindow window = new IsolationWindow(peptide, targetMz, charge, rtMin, rtMax, isDecoy);
				isolationWindows.add(window);
				
				if (!hasPrintedAddingPrecursorInfo) {
					hasPrintedAddingPrecursorInfo = true;
					System.out.println("Adding an mz at " + targetMz + " and RT " + rtCenter + " min " + rtMin/60 + " max " + rtMax/60 );

				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}
	
	// This function is meant to detect the delimiter - Thermo Mass lists are usually in a .csv when exported from the method editor, but TargetedBootstrapper provides .txt. This will let us accept either for now. 
	private static String getDelimiter(String filePath) {
		String lowerPath = filePath.toLowerCase();
		
		if (lowerPath.endsWith(".csv")) {
			return ",";
		}
		
		if (lowerPath.endsWith(".txt")) {
			return "\t";
		}
		
		if (lowerPath.endsWith(".tsv")) {
			return "\t";
		}
		
		throw new IllegalArgumentException("Error in reading IsolationWindows. Mass List file must be .csv, .txt, or .tsv. Current file is: " + filePath);

	}
}
