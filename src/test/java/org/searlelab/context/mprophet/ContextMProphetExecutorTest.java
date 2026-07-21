package org.searlelab.context.mprophet;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;
import org.searlelab.context.mprophet.ContextMProphetExecutor;


public class ContextMProphetExecutorTest {

	@Test
	
	public void test() {

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

}
