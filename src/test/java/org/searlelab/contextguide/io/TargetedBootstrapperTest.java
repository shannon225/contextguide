package org.searlelab.contextguide.io;

import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.Objects;

import org.junit.Test;
import org.junit.Rule;
import org.junit.rules.TemporaryFolder;
import org.searlelab.context.io.TargetedBootstrapper;

public class TargetedBootstrapperTest {

	@Rule
	public TemporaryFolder tempFolder = TemporaryFolder.builder().assureDeletion().build();

	@Test
	public void testIfBootstrapperCreatesExpectedFiles() throws Throwable {
		Path testDirectory = tempFolder.newFolder("targeted-bootstrapper-test").toPath();

		URL libraryURL = Objects.requireNonNull(getClass().getClassLoader().getResource("IL2_and_IL15_Combo.elib"));
		String library = Paths.get(libraryURL.toURI()).toString();

		URL diaURL = Objects.requireNonNull(getClass().getClassLoader().getResource("IL2A_GPFDIA_0combined_masked0_assay.dia"));
		Path sourceDIAPath = Paths.get(diaURL.toURI());

		Path testDIAPath = testDirectory.resolve("IL2A_GPFDIA_0combined.dia");

		Files.copy(sourceDIAPath, testDIAPath, StandardCopyOption.REPLACE_EXISTING);
		Path targetDecoyMap = testDirectory.resolve("target_decoy_map.txt");
		TargetedBootstrapper bootstrapper = new TargetedBootstrapper();
		bootstrapper.execute(library.toString(), testDIAPath.toString(), targetDecoyMap, 0, 10, 0.5f, 1.0);

	}
}