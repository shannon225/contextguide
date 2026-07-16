package org.searlelab.contextguide.mprophet;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.searlelab.context.mprophet.ContextFeatureScorer;
import org.searlelab.context.mprophet.ScoredFeature;

public class ContextFeatureScorerTest {

	@Rule
	public TemporaryFolder tempFolder = TemporaryFolder.builder().assureDeletion().build();
	
	@Test
	public void test() throws Exception {
		// This tests if the ContextFeatureScorer can find the correct files and return scored features 
		
		// Locate the files for the test	
		URL rawFileName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.dia");
		URL libraryFileName = getClass().getClassLoader().getResource("IL2_and_IL15_Combo.elib");
		URL fastaFileName = getClass().getClassLoader().getResource("mus_musculus_reviewed_uniprot.fasta");
		URL massListFileName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.txt");

		// These files must exist for the test to pass
		assertNotNull("DIA file was not found. ", rawFileName);
		assertNotNull("Library was not found. ", libraryFileName);
		assertNotNull("Fasta was not found. ", fastaFileName);
		assertNotNull("Mass list was not found. ", massListFileName);
		
		// Copy the files to the temporary directory

		Path rawFilePath = Paths.get(rawFileName.toURI());
		Path libraryFilePath = Paths.get(libraryFileName.toURI());
		Path fastaFilePath = Paths.get(fastaFileName.toURI());
		
		String massListPath = Paths.get(massListFileName.toURI()).toString();

		File rawFile = rawFilePath.toFile(); 
		File library = libraryFilePath.toFile();
		File fasta = fastaFilePath.toFile();
		
		String baseName = rawFile.toString().replaceFirst("\\.dia$", "");

			ArrayList<ScoredFeature> partitionedFeatures = ContextFeatureScorer.scoreFeatures(library, rawFile, fasta, baseName, massListPath);
			assertNotNull(partitionedFeatures);
			
			System.out.println("Last feature in the list: " + partitionedFeatures.getLast() + " peptides remaining.");
			
			// Verify that the expected output files were created
		Path referenceOutput = Paths.get(baseName + "_reference.features.txt");
		Path backgroundOutput = Paths.get(baseName + "_background.features.txt");
		
		assertTrue("Reference feature file was not created.", Files.exists(referenceOutput));
		assertTrue("Background feature file was not created.", Files.exists(backgroundOutput));
		
		System.out.println(partitionedFeatures.size() + " paritioned features were returned.");
		int featuresSize = 17;
		assertTrue(partitionedFeatures.size()==featuresSize);
		

		}
	}


