package org.searlelab.contextguide.mprophet;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Paths;
import java.util.ArrayList;
import org.junit.Test;
import org.searlelab.contextguide.mprophet.ContextFeatureScorer;

public class ContextFeatureScorerTest {
// TODO This test is too long. It processes a whole file. Make this test run more quickly by reducing of features in the file to 10 
	@Test
	public void smokeTest() throws Exception {
		
		URL rawFileName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_10kb.dia");
		URL libraryFileName = getClass().getClassLoader().getResource("IL2_and_IL15_Combo.elib");
		URL fastaFileName = getClass().getClassLoader().getResource("mus_musculus_reviewed_uniprot.fasta");
		URL massListFileName = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.txt");

		String rawFilePath = Paths.get(rawFileName.toURI()).toString();
		String libraryFilePath = Paths.get(libraryFileName.toURI()).toString();
		String fastaFilePath = Paths.get(fastaFileName.toURI()).toString();
		String massListPath = Paths.get(massListFileName.toURI()).toString();

		
		
//		String rawFilePath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_5min/IL2A_GPFDIA_0combined_Masked0_assay.dia";
//		String libraryFilePath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_5min/IL2_and_IL15_Combo.elib";
//		String fastaPath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_5min/mus_musculus_reviewed_uniprot.fasta";
//		String massListPath = "C:/Users/m334793/Documents/asms2026/stellar/2mz_5min/IL2A_GPFDIA_0combined_masked0_assay.txt";

		String baseName = rawFilePath.replaceFirst("\\.dia$", "");
		final File fasta = new File(fastaFilePath);
		File rawFile = new File(rawFilePath);
		File library = new File(libraryFilePath);

		try {			
			ArrayList<ScoredFeature> partitionedFeatures = ContextFeatureScorer.scoreFeatures(library, rawFile, fasta, baseName, massListPath);
			assertNotNull(partitionedFeatures);
			
			System.out.println("Last feature in the list: " + partitionedFeatures.getLast());
			
		//	ScoredFeature testFeature = ScoredFeature();
	//		assertEquals(partitionedFeatures.getLast(), );
		} catch (Exception e) {
			System.out.println("Something did not work when trying to score features... see the error tace");
			e.printStackTrace();
		} finally {
		System.out.println("Features have been scored.");
		
		}
		

		}
	}


