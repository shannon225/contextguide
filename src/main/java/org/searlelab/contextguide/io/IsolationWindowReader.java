package org.searlelab.contextguide.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.contextguide.mprophet.IsolationWindow;

public class IsolationWindowReader {

//	private static final String DELIM = "," || "\t"; 
//	public static void main(String[] args) {

	//	String massListFile = "C:/Users/m334793/Documents/Library/assay7.csv";

		//ArrayList<IsolationWindow> isolationWindows = parseMassList(massListFile);
		
		//System.out.println("Total windows read: " + isolationWindows.size());

	//}

	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
			String delim = getDelimiter(massListFile);
	
			@SuppressWarnings("unused")
			String header = br.readLine();
			
			String line;
			while ((line = br.readLine()) != null) {

				String[] columns = line.split(delim, -1);
				
				 boolean hasPrintedDebugInfo = false;
				 if (!hasPrintedDebugInfo) {
				        System.out.println("Line being read: " + line);
				        System.out.println("Number of columns: " + columns.length);
				        hasPrintedDebugInfo = true;
				    }

			//	String columns[] = line.split(DELIM, -1);
				System.out.println(line); // Console will print what the data looks like as its read in

				String peptide = columns[0];
				// adduct = 1
				// compound = 2
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
				System.out.println("Adding an mz at " + targetMz + " and RT " + rtCenter + " min " + rtMin/60 + " max " + rtMax/60 
//						+ "\nRTCenter: " + rtCenter
//						+ "\ntargetMz: " + targetMz
//						+ "\nrtStart: " + rtMin
//						+ "\nrtStop: " + rtMax
						);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}
	
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

	    throw new IllegalArgumentException(
	        "Mass list file must be a .csv or .tsv file: " + filePath
	    );
	}
}
