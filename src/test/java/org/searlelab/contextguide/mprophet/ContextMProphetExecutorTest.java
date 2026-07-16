package org.searlelab.contextguide.mprophet;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.searlelab.context.mprophet.ContextMProphetExecutor;

import com.google.common.io.Files;

public class ContextMProphetExecutorTest {

	@Rule 
	public TemporaryFolder tempFolder = TemporaryFolder.builder().assureDeletion().build();
	
	
	
	@Test
	public void test() throws Throwable {

		// Map files 
		URL libraryFileName = getClass().getClassLoader().getResource("IL2_and_IL15_Combo.elib");
		URL fastaFileName = getClass().getClassLoader().getResource("mus_musculus_reviewed_uniprot.fasta");

		// Where the feature files are located: 
		URL diaFileName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.dia");
		URL massListName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.txt");

		// Are all the paths mapped to real files? 
		assertNotNull("Library was not found. ", libraryFileName);
		assertNotNull("Fasta was not found. ", fastaFileName);
		assertNotNull("DIA file was not found. ", diaFileName);
		assertNotNull("Mass list was not found. ", massListName);

		Path libraryPath = Paths.get(libraryFileName.toURI());
		Path fastaPath = Paths.get(fastaFileName.toURI());
		Path diaPath = Paths.get(diaFileName.toURI());
		Path massListPath = Paths.get(massListName.toURI());
		
		File copiedLibraryFile =   new File(tempFolder.getRoot(), libraryPath.getFileName().toString());
		Files.copy(libraryPath.toFile(), copiedLibraryFile);
		
		File diaFile = diaPath.toFile();
		//		File[] diaFiles = diaFolder.listFiles();

		//		System.out.println("DIA Folder was identified as: " + diaFolder.getAbsolutePath());
		//		System.out.println("DIA files detected! The following files will be processed with MProphet:" + diaFolder.listFiles());


		// Ignore files that do not end in .dia 
		if (!diaFile.isFile() && !diaFile.getName().endsWith(".dia")) {
			
		String diaName = diaFile.getName(); 
		String baseName = diaName.substring(0, diaName.lastIndexOf(".dia"));

		File massListFile = new File(diaFile, baseName + ".txt");

		System.out.println("Processessing " + diaFile.getName());

		if (!massListFile.exists()) {
			System.out.println("Skipping " + diaFile.getName() + " because mass list was not found: " + massListPath);
			
		}
		ContextMProphetExecutor.executeContextMProphet(libraryPath.toString(), fastaPath.toString(), diaPath.toString(), massListPath.toString());
	}	
}
	
}

