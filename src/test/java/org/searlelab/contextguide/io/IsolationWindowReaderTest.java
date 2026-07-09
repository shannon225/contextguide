package org.searlelab.contextguide.io;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.searlelab.contextguide.io.IsolationWindowReader;
import org.searlelab.contextguide.mprophet.IsolationWindow;
import org.junit.Test;

public class IsolationWindowReaderTest {

	@Test
	public void testThatFileWorks() throws Throwable {
	// Does this correctly read a mass list? 
		URL massListFile = getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.txt");
			
		String massListPath = Paths.get(massListFile.toURI()).toString();
		System.out.println("Mass list test file is " + massListPath);
		ArrayList<IsolationWindow> massList = IsolationWindowReader.parseMassList(massListPath);
		
		assertNotNull(massList);
	}

}
