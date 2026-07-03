package org.searlelab.contextguide.mprophet;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;

public class ContextWindowExtractor {
	
	public static void main(String[] args) {
		Path targetAssayPath = Paths.get("C:/Users/m334793/git/Documents/Library/assay7.csv");
		Path rawFilePath = Paths.get("C/Users/m334794/git/Documents/Library/masked1_cd14_combined.dia");
		Path outputPath = Paths.get("C:/Users/m334793/git/Documents/Library");
	
//		if (!Files.exists(targetAssayPath)) {
//			System.out.println("Target assay does not exist.");
//		}
		
//		if (!Files.isRegularFile(targetAssayPath)) {
//			System.out.println("Target assay path is not a file. Ensure it is a .dia/.mzML/.raw file.");
//		}
		
//		if (!Files.isReadable(targetAssayPath)) {
//			System.out.println("Target assay could not be read.");
//		}
	}
	
	public static StripeFile maskByTargetWindows(File rawFilePath,  String outputPath, IsolationWindow isolationWindows) throws IOException {
		long startTime = System.nanoTime();
		
		// Open the masked File
		StripeFile maskedFile = new StripeFile(false);
		
		// Make a file output from the argument 'output path' 
		File outputFile = new File(outputPath);
		StripeFile rawFile = new StripeFile(false);

		// Read in isolation windows list
		double targetMz;
		float rtMin;
		float rtMax;
		boolean isDecoy;
		

		try {
			rawFile.openFile(rawFilePath);
			
			// Add ranges
			HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();
			System.out.println("Masking DIA file based on assay mass list.");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		long stopTime = System.nanoTime();
		long duration = stopTime - startTime;
		System.out.println("maskByTargetWindows(): Time taken (ms) : " + duration / 1_000_000);

		return maskedFile;

	}
	
//	public static maskByTargetMz() {
		
//	}

}
