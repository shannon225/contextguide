package org.searlelab.context.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.context.mprophet.IsolationWindow;

public class IsolationWindowReader {

	private static final String DELIM = ","; // How do I make this a csv?

	public static void main(String[] args) {

		String massListFile = "C:/Users/m334793/Documents/Library/assay7.csv";

		ArrayList<IsolationWindow> isolationWindows = parseMassList(massListFile);
		
		System.out.println("Total windows read: " + isolationWindows.size());

	}

	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
	
			@SuppressWarnings("unused")
			String header = br.readLine();
			
			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split(DELIM, -1);
				System.out.println("Added another line " + line); // Console will print what the data looks like as its read in

				double targetMz = Double.parseDouble(columns[3]);
				float rtCenter = Float.parseFloat(columns[5]);
				float rtWindow = Float.parseFloat(columns[6]);

				float rtMin = (rtCenter - (rtWindow / 2))*60;
				float rtMax = (rtCenter + (rtWindow / 2))*60;

				boolean isDecoy = false;

				// Assemble each window
				IsolationWindow window = new IsolationWindow(targetMz, rtMin, rtMax, isDecoy);
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
}
