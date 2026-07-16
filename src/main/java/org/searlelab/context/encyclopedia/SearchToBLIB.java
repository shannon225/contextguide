package org.searlelab.context.encyclopedia;

import static edu.washington.gs.maccoss.encyclopedia.Encyclopedia.QUIET_MODE_ARG;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.DataFormatException;

import org.searlelab.msrawjava.io.StripeFileInterface;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.primitives.Floats;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.ScribeTwo;
import edu.washington.gs.maccoss.encyclopedia.SearchGUIMain;
import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper.PickedProteinFDRResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AbstractRetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.EncyclopediaTwoPeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface.AlignmentDataPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.SimplePeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.JavaPotExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetReiter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.LocalizationDataToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.precursor.DDAPrecursorIntegrator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.IntensityNormalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeTwoDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantXCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DDASearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.OverlappingDiaPreprocessor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.ParsingUtils;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.FileLogRecorder;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.ThrowingConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.iterator.TDoubleIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.hash.TDoubleHashSet;

public class SearchToBLIB {
	public static final String GLOBAL_NAME = "global";

	public static void main(String[] args) {
		final Pair<List<String>, HashMap<String, String>> parsedArgs = CommandLineParser.parseMultipleAndRemainingArguments(args, Encyclopedia.INPUT_DIA_TAG);
		final List<String> diaPaths = parsedArgs.x;
		HashMap<String, String> arguments = parsedArgs.y;
		arguments=InstrumentSpecificSearchParameters.checkParameters(arguments);

		if (arguments.size()==0) {
			SearchGUIMain.runGUI(ProgramType.EncyclopeDIA);
		} else if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("SearchToLIB Help");
			Logger.timelessLogLine("You should prefix your arguments with a high memory setting, e.g. \"-Xmx8g\" for 8gb");
			Logger.timelessLogLine("Other Programs: ");
			Logger.timelessLogLine("\t-pecan\trun Pecanpie export (use -pecan -h for Pecan help)");
			Logger.timelessLogLine("\t-xcordia\trun XCorDIA export (use -xcordia -h for XCorDIA help)");
			Logger.timelessLogLine("\t-phospho\trun phospho localization export (use -phospho -h for localization help)");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .DIA or .MZML file or directory");
			Logger.timelessLogLine("\t-o\toutput library .ELIB file");
			Logger.timelessLogLine("\t-a\talign between files (default=true)");
			Logger.timelessLogLine("\t-blib\twrite .BLIB instead of .ELIB (default=false)");
			Logger.timelessLogLine("Potentially Required Parameters: ");
			Logger.timelessLogLine("\t-l\toriginal searched library .DLIB or .ELIB file (required by EncyclopeDIA Export)");
			Logger.timelessLogLine("\t-f\toriginal fasta file (required by Pecan/XCorDIA Export)");
			Logger.timelessLogLine("\t-t\toriginal target file (optional for Pecan/XCorDIA Export)");

			Logger.timelessLogLine("Other Parameters: ");
			TreeMap<String, String> defaults=new TreeMap<String, String>(SearchParameterParser.getExportParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			System.exit(1);
			
		} else if (arguments.containsKey("-v")||arguments.containsKey("-version")||arguments.containsKey("--version")) {
			Logger.logLine("EncyclopeDIA SearchToLIB version "+ProgramType.getGlobalVersion().toString());
			System.exit(1);
			
		} else {
			if (diaPaths.isEmpty()) {
				Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
				System.exit(1);
			}

			final List<File> diaFiles = diaPaths.stream()
					.map(Paths::get)
					.map(Path::toFile)
					.collect(Collectors.toList());

			try {
				if (arguments.containsKey("-pecan")||arguments.containsKey("-walnut")) {
					VersioningDetector.checkVersionCLI(ProgramType.PecanPie);
					convertPecan(diaFiles, arguments);
				} else if (arguments.containsKey("-xcordia")) {
					VersioningDetector.checkVersionCLI(ProgramType.XCorDIA);
					convertXCorDIA(diaFiles, arguments);
				} else if (arguments.containsKey("-scribe")) {
					convertScribe(diaFiles, arguments);

				} else {
					VersioningDetector.checkVersionCLI(ProgramType.EncyclopeDIA);
					convertEncyclopedia(diaFiles, arguments);
				}

			} catch (Exception e) {
				Logger.errorLine("Encountered Fatal Error!");
				Logger.errorException(e);

				// Forcibly exit with error to avoid hanging the process if there are any leftover user threads
				// that weren't properly cleaned up as a result of the error. First we log any hung / remaining
				// threads though, to aid debugging.

				for (Entry<Thread, StackTraceElement[]> threadEntry : Thread.getAllStackTraces().entrySet()) {
					if (threadEntry.getKey().isDaemon()) {
						continue;
					}
					if (Thread.currentThread().equals(threadEntry.getKey())) {
						continue;
					}

					Logger.errorLine("\nLeftover user thread will be KILLED: " + threadEntry.getKey().getName());
					for (StackTraceElement element : threadEntry.getValue()) {
						Logger.errorLine("  " + element);
					}
				}

				System.exit(1); // nonzero status indicates error
			} finally {
				Logger.close();
			}
		}
	}

	public static void convertScribe(List<File> diaFiles, HashMap<String, String> arguments) {
		if (!arguments.containsKey("-l")||!arguments.containsKey("-o")||!arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input library file (-l), a fasta database (-f), and an output library file (-o)");
			System.exit(1);
		}

		arguments=InstrumentSpecificSearchParameters.checkParameters(arguments);

		File fastaFile=new File(arguments.get("-f"));
		File libraryFile=new File(arguments.get("-l"));
		File outputFile=new File(arguments.get("-o"));

		final boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		final boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);
		final boolean alignOnly = ParsingUtils.getBoolean("-alignOnly", arguments, false);

		final SearchParameters parameters=SearchParameterParser.parseParameters(arguments);

		final OutputFormat outputFormat;

		if (!alignOnly) {
			outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;
		} else {
			if (!alignBetweenFiles) {
				Logger.errorLine("-alignOnly requires alignment to be enabled; try running with `-a true`");
				System.exit(1);
			}

			if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
				Logger.errorLine("-alignOnly requires -quantifyAcrossSamples true");
				System.exit(1);
			}

			if (writeBlib) {
				Logger.errorLine("-alignOnly requires ELIB output; try running with `-blib false`");
				System.exit(1);
			}

			if (arguments.containsKey("-alignmentFrom")) {
				Logger.errorLine("Error: -alignOnly and -alignmentFrom are incompatible");
				System.exit(1);
			}

			outputFormat = OutputFormat.ALIB;
		}

		ScribeScoringFactory factory = new ScribeScoringFactory(parameters);
		Logger.timelessLogLine("SearchToLIB EncyclopeDIA version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		for (File diaFile : diaFiles) {
			Logger.timelessLogLine(" -i " + diaFile.getAbsolutePath());
		}
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -l "+libraryFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(" -alignOnly " + alignOnly);
		Logger.timelessLogLine(parameters.toString());

		if (arguments.containsKey(QUIET_MODE_ARG)) {
			Logger.PRINT_TO_SCREEN = false;
		}

		try {
			LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			for (File diaFile: diaFiles) {
				if (diaFile.isDirectory()) {
					File[] files = diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
					if (files.length == 0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files: " + diaFile.getAbsolutePath());
						System.exit(1);
					}

					for (File file : files) {
						ScribeJobData job = new ScribeJobData(file, fastaFile, library, factory);
						pecanJobs.add(job);
					}
				} else if (alignOnly && !diaFile.exists()) {
					// Special case -- when running alignment-only we may not have the .DIA available but want
					// to handle the job using Percolator/ELIB results only.
					//pecanJobs.add(EncyclopediaJobData.getDummyFor(diaFile, fastaFile, library, factory));
					
					// FIXME: this edge case does not work! Throw error instead
					Logger.errorLine("Unexpected mode running Scribe quantification with alignment-only data without .DIA files available! DIA files are required for Scribe quant.");
					System.exit(1);
				} else {
					ScribeJobData job = new ScribeJobData(diaFile, fastaFile, library, factory);
					pecanJobs.add(job);
				}
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");

			if (!arguments.containsKey("-alignmentFrom")) {
				// Main program: convert to appropriate format
				convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles, parameters);
			} else {
				// Sub-program: quantify from previously-computed alignment/transition refinement

				if (!alignBetweenFiles) {
					Logger.errorLine("-alignmentFrom requires alignment to be enabled; try running with `-a true`");
					System.exit(1);
				}

				if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
					Logger.errorLine("-alignmentFrom requires -quantifyAcrossSamples true");
					System.exit(1);
				}

				if (writeBlib) {
					Logger.errorLine("-alignmentFrom requires ELIB output; try running with `-blib false`");
					System.exit(1);
				}

				if (alignOnly) {
					Logger.errorLine("Error: -alignOnly and -alignmentFrom are incompatible");
					System.exit(1);
				}

				convertElibQuantOnly(new EmptyProgressIndicator(), pecanJobs, outputFile, new File(arguments.get("-alignmentFrom")), parameters);
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convertXCorDIA(List<File> diaFiles, HashMap<String, String> arguments) {
		if (!arguments.containsKey("-f")||!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input fasta file (-f) and an output library file (-o)");
			System.exit(1);
		}

		File fastaFile=new File(arguments.get("-f"));
		File outputFile=new File(arguments.get("-o"));
		boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);

		final OutputFormat outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;

		PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
		XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(parameters);
		Logger.timelessLogLine("SearchToLIB XCorDIA version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		for (File diaFile : diaFiles) {
			Logger.timelessLogLine(" -i " + diaFile.getAbsolutePath());
		}
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(parameters.toString());

		try {
			ArrayList<FastaPeptideEntry> targets;
			if (arguments.containsKey(XCorDIA.TARGET_FASTA_TAG)) {
				targets=FastaReader.readPeptideFasta(new File(arguments.get(XCorDIA.TARGET_FASTA_TAG)), parameters);
			} else {
				targets=null;
			}
			LibraryInterface library;
			if (arguments.containsKey("-l")) {
				library=BlibToLibraryConverter.getFile(new File(arguments.get("-l")));
			} else {
				library=null;
			}
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			for (File diaFile : diaFiles) {
				if (diaFile.isDirectory()) {
					File[] files = diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
					if (files.length == 0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files: " + diaFile.getAbsolutePath());
						System.exit(1);
					}

					for (File file : files) {
						XCorDIAJobData job = new XCorDIAJobData(Optional.ofNullable(targets), Optional.ofNullable(library), file, fastaFile, factory);
						pecanJobs.add(job);
					}
				} else {
					XCorDIAJobData job = new XCorDIAJobData(Optional.ofNullable(targets), Optional.ofNullable(library), diaFile, fastaFile, factory);
					pecanJobs.add(job);
				}
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");
			convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles, parameters);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convertPecan(List<File> diaFiles, HashMap<String, String> arguments) {
		if (!arguments.containsKey("-f")||!arguments.containsKey("-o")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input fasta file (-f) and an output library file (-o)");
			System.exit(1);
		}

		File fastaFile=new File(arguments.get("-f"));
		File outputFile=new File(arguments.get("-o"));
		boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);

		final OutputFormat outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;

		PecanSearchParameters parameters=PecanParameterParser.parseParameters(arguments);
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, outputFile);
		Logger.logLine("SearchToLIB Pecan version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		for (File diaFile : diaFiles) {
			Logger.timelessLogLine(" -i " + diaFile.getAbsolutePath());
		}
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(parameters.toString());

		try {
			ArrayList<FastaPeptideEntry> targets;
			if (arguments.containsKey(Pecanpie.TARGET_FASTA_TAG)) {
				targets=FastaReader.readPeptideFasta(new File(arguments.get(Pecanpie.TARGET_FASTA_TAG)), parameters);
			} else {
				targets=null;
			}
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			for (File diaFile : diaFiles) {
				if (diaFile.isDirectory()) {
					File[] files = diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
					if (files.length == 0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files: " + diaFile.getAbsolutePath());
						System.exit(1);
					}

					for (File file : files) {
						PecanJobData job = new PecanJobData(Optional.ofNullable(targets), file, fastaFile, factory);
						pecanJobs.add(job);
					}
				} else {
					PecanJobData job = new PecanJobData(Optional.ofNullable(targets), diaFile, fastaFile, factory);
					pecanJobs.add(job);
				}
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");
			convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles, parameters);
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public static void convertEncyclopedia(List<File> diaFiles, HashMap<String, String> arguments) {
		if (!arguments.containsKey("-l")||!arguments.containsKey("-o")||!arguments.containsKey("-f")) {
			Logger.errorLine("You are required to specify an input file or directory (-i), an input library file (-l), a fasta database (-f), and an output library file (-o)");
			System.exit(1);
		}
		
		arguments=InstrumentSpecificSearchParameters.checkParameters(arguments);

		File fastaFile=new File(arguments.get("-f"));
		File libraryFile=new File(arguments.get("-l"));
		File outputFile=new File(arguments.get("-o"));

		final boolean alignBetweenFiles=ParsingUtils.getBoolean("-a", arguments, true);
		final boolean writeBlib=ParsingUtils.getBoolean("-blib", arguments, false);
		final boolean alignOnly = ParsingUtils.getBoolean("-alignOnly", arguments, false);

		final SearchParameters parameters=SearchParameterParser.parseParameters(arguments);

		final OutputFormat outputFormat;

		if (!alignOnly) {
			outputFormat = writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB;
		} else {
			if (!alignBetweenFiles) {
				Logger.errorLine("-alignOnly requires alignment to be enabled; try running with `-a true`");
				System.exit(1);
			}

			if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
				Logger.errorLine("-alignOnly requires -quantifyAcrossSamples true");
				System.exit(1);
			}

			if (writeBlib) {
				Logger.errorLine("-alignOnly requires ELIB output; try running with `-blib false`");
				System.exit(1);
			}

			if (arguments.containsKey("-alignmentFrom")) {
				Logger.errorLine("Error: -alignOnly and -alignmentFrom are incompatible");
				System.exit(1);
			}

			outputFormat = OutputFormat.ALIB;
		}

		LibraryScoringFactory factory=EncyclopediaScoringFactory.getScoringFactory(arguments, parameters);
		Logger.timelessLogLine("SearchToLIB EncyclopeDIA version "+ProgramType.getGlobalVersion().toString());

		Logger.timelessLogLine("Parameters:");
		for (File diaFile : diaFiles) {
			Logger.timelessLogLine(" -i " + diaFile.getAbsolutePath());
		}
		Logger.timelessLogLine(" -f "+fastaFile.getAbsolutePath());
		Logger.timelessLogLine(" -l "+libraryFile.getAbsolutePath());
		Logger.timelessLogLine(" -o "+outputFile.getAbsolutePath());
		Logger.timelessLogLine(" -a "+alignBetweenFiles);
		Logger.timelessLogLine(" -blib "+writeBlib);
		Logger.timelessLogLine(" -alignOnly " + alignOnly);
		Logger.timelessLogLine(parameters.toString());

		if (arguments.containsKey(QUIET_MODE_ARG)) {
			Logger.PRINT_TO_SCREEN = false;
		}

		try {
			LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
			
			ArrayList<SearchJobData> pecanJobs=new ArrayList<SearchJobData>();
			for (File diaFile: diaFiles) {
				if (diaFile.isDirectory()) {
					File[] files = diaFile.listFiles(StripeFileGenerator.getFilenameFilter());
					if (files.length == 0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files: " + diaFile.getAbsolutePath());
						System.exit(1);
					}

					if (files.length == 0) {
						Logger.errorLine("Your specified input (-i) directory didn't contain any .RAW files!");
						System.exit(1);
					}
					for (File file : files) {
						EncyclopediaJobData job = new EncyclopediaJobData(file, fastaFile, library, factory);
						pecanJobs.add(job);
					}
				} else if (alignOnly && !diaFile.exists()) {
					// Special case -- when running alignment-only we may not have the .DIA available but want
					// to handle the job using Percolator/ELIB results only.
					pecanJobs.add(EncyclopediaJobData.getDummyFor(diaFile, fastaFile, library, factory));
				} else {
					EncyclopediaJobData job = new EncyclopediaJobData(diaFile, fastaFile, library, factory);
					pecanJobs.add(job);
				}
			}
			Logger.logLine("Attempting to process "+pecanJobs.size()+" searches...");

			if (!arguments.containsKey("-alignmentFrom")) {
				// Main program: convert to appropriate format
				convert(new EmptyProgressIndicator(), pecanJobs, outputFile, outputFormat, alignBetweenFiles, parameters);
			} else {
				// Sub-program: quantify from previously-computed alignment/transition refinement

				if (!alignBetweenFiles) {
					Logger.errorLine("-alignmentFrom requires alignment to be enabled; try running with `-a true`");
					System.exit(1);
				}

				if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
					Logger.errorLine("-alignmentFrom requires -quantifyAcrossSamples true");
					System.exit(1);
				}

				if (writeBlib) {
					Logger.errorLine("-alignmentFrom requires ELIB output; try running with `-blib false`");
					System.exit(1);
				}

				if (alignOnly) {
					Logger.errorLine("Error: -alignOnly and -alignmentFrom are incompatible");
					System.exit(1);
				}

				convertElibQuantOnly(new EmptyProgressIndicator(), pecanJobs, outputFile, new File(arguments.get("-alignmentFrom")), parameters);
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}

	public enum OutputFormat {
		/**
		 * Write to the ELIB format. If {@code inferrer} is present, the resulting file will be a "quantitative" ELIB,
		 * using the precomputed top-N transitions for quantification and inferred (aligned) RTs when the peptide was
		 * not detected in the initial single-file search. Additionally, quantitative matrices for peptides and proteins
		 * will be written.
		 *
		 * All passing peptides will be included.
		 *
		 * {@code globalPercolatorPeptides} should be provided if converting more than a single search, but should otherwise be empty.
		 */
		ELIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters, boolean integratePrecursors) {
				convertElib(progress, jobs, outputFile, Optional.of(passingPeptides), globalPercolatorFiles, inferrer, parameters, integratePrecursors);
			}
		},

		/**
		 * Write results to the "alignment-only library" format, which records the passing peptides, RT alignment, and
		 * refined transitions for the experiment to a library file without performing any additional work. The resulting
		 * file can then be used to quantify the same targets in later separate runs of one (or more) sample(s).
		 *
		 * Note that {@code inferrer} must be present to support this export type.
		 */
		ALIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters, boolean integratePrecursors) {
				if (!inferrer.isPresent()) {
					throw new IllegalArgumentException("Unable to export alignment-only library without RT alignment and transition refinement!");
				}

				convertAlib(progress, jobs, outputFile, passingPeptides, globalPercolatorFiles, inferrer.get(), parameters);
			}
		},

		/**
		 * Write to the BLIB format, suitable for use with Skyline. Only quantifiable peptides will be written from each
		 * search. Additionally, a TSV "integration" file will be written with details of the peptides included in the library.
		 */
		BLIB {
			@Override
			void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters, boolean integratePrecursors) {
				convertBlib(progress, jobs, outputFile, Optional.of(passingPeptides.x), inferrer);
			}
		};

		/**
		 * Write results to the given location in this format. Typically this method should only be called from
		 * {@link SearchToBLIB#convert(ProgressIndicator, List, File, OutputFormat, boolean, SearchParameters)}
		 * which will handle either reading or computing the necessary information for a group of samples.
		 *
		 * Will also compute and output related information in some cases, depending on the format.
		 *
		 * @param progress A progress indicator that will be used during the conversion process
		 * @param jobs The jobs whose results should be included in the output file
		 * @param outputFile The location where the new library will be created (will be overwritten if it exists)
		 * @param passingPeptides The results of running Percolator to determine the list of peptides that will be
		 *                        included, as returned by {@link PercolatorReader#getPassingPeptidesFromTSV}
		 * @param globalPercolatorFiles Used by some formats to get additional information when Percolator has been run on
		 *                              results from multiple input files
		 * @param inferrer If aligning between files, the inferrer which provides RT alignment and consistent, refined transitions
		 * @param parameters The parameters that should be used during conversion and (in some cases) written to the output file
		 */
		abstract void convert(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters, boolean integratePrecursors);
	}

	/**
	 * Legacy form of {@link #convert(ProgressIndicator, List, File, OutputFormat, boolean, SearchParameters)}
	 * which supports only ELIB and BLIB formats.
	 *
	 * @see #convert(ProgressIndicator, List, File, OutputFormat, boolean, SearchParameters)
	 *
	 * @deprecated it's better to directly specify the desired output format with an enum constant
	 */
	@Deprecated
	public static void convert(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File libFile, boolean writeBlib, boolean alignBetweenFiles) {
		convert(
				progress,
				pecanJobs,
				libFile,
				writeBlib ? OutputFormat.BLIB : OutputFormat.ELIB,
				alignBetweenFiles,
				null
		);
	}

	/**
	 * For the given previously-run single-file searches (jobs), gather or compute the necessary information to create
	 * a combined output in the given format. This handles the core jobs of (if necessary) running Percolator, reading
	 * Percolator results to determine the set of global passing peptides, performing retention time alignment and
	 * transition refinement (if {@code alignBetweenFiles} is true), and writing results to the given output file, which
	 * may involve additional work like quantifying peptides in each single file.
	 *
	 * @param progress A progress indicator that will be used during the conversion process
	 * @param origPecanJobs The jobs whose results should be included in the output file
	 * @param libFile The location where the new library will be created (will be overwritten if it exists)
	 * @param outputFormat The format which should be written
	 * @param alignBetweenFiles If RT alignment
	 */
	public static void convert(ProgressIndicator progress, List<? extends SearchJobData> origPecanJobs, File libFile, OutputFormat outputFormat, boolean alignBetweenFiles, SearchParameters parameters) {
		ArrayList<SearchJobData> preprocessedPecanJobs=OverlappingDiaPreprocessor.preprocess(progress, new ArrayList<SearchJobData>(origPecanJobs));
		List<? extends SearchJobData> pecanJobs = getProcessedJobs(preprocessedPecanJobs);
		if (pecanJobs.size()==0) {
			Logger.errorLine("Can't find any representative jobs! Failing...");

			for (int i=0; i<origPecanJobs.size(); i++) {
				SearchJobData job=origPecanJobs.get(i);
				Logger.errorLine(" Checking raw file "+(i+1)+": "+job.getDiaFileReader().getFile().exists());
				Logger.errorLine(" Checking feature file "+(i+1)+": "+job.getPercolatorFiles().getInputTSV().exists());
				Logger.errorLine(" Checking result file "+(i+1)+": "+job.getPercolatorFiles().getPeptideOutputFile().exists());
			}
			return;
		}

		final SearchJobData representativeJob=pecanJobs.get(0);
		
		ArrayList<File> featureFiles=new ArrayList<File>();
		for (int i=0; i<pecanJobs.size(); i++) {
			SearchJobData job=pecanJobs.get(i);
			featureFiles.add(job.getPercolatorFiles().getInputTSV());
		}

		Logger.logLine("Using "+representativeJob.getOriginalDiaFileName()+" to extract representative search parameters");
		if (parameters==null) {
			Logger.logLine("Using "+representativeJob.getOriginalDiaFileName()+" to extract representative search parameters");
			parameters=representativeJob.getParameters();
		}
		
		boolean integratePrecursors=parameters.isIntegratePrecursors();

		boolean useTDC=false;
		boolean anyDDA=false;
		for (int i=0; i<pecanJobs.size(); i++) {
			SearchJobData job=pecanJobs.get(i);
			if (job instanceof DDASearchJobData) {
				anyDDA=true;
				useTDC=true;
				break;
			}
			if (job instanceof EncyclopediaTwoJobData||job instanceof ScribeTwoDIAJobData) {
				useTDC=true;
				break;
			}
		}
		if (useTDC) {
			Logger.logLine("Running Percolator in target-decoy completition mode");
		}
		if (anyDDA) {
			if (!integratePrecursors) {
				Logger.logLine("Found DDA data, forcing integration of precursors");
				integratePrecursors=true;
			}
		} else {
			Logger.logLine("Found DIA data only, running integration of fragments");
		}

		String filename=libFile.getName();
		if (filename.lastIndexOf('.')>0) {
			filename=filename.substring(0, filename.lastIndexOf('.'));
		}
		File bigFeatureFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_results.txt");
		File bigPercolatorDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_decoy.txt");
		File bigPercolatorProteinFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_results.txt");
		File bigPercolatorProteinDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_decoy.txt");
		PercolatorExecutionData bigPercolatorFiles=new PercolatorExecutionData(bigFeatureFile, representativeJob.getPercolatorFiles().getFastaFile(), 
				bigPercolatorFile, bigPercolatorDecoyFile, bigPercolatorProteinFile, bigPercolatorProteinDecoyFile, parameters, !useTDC);

		final float threshold=parameters.getEffectivePercolatorThreshold();
		try {
			Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Boolean> percolatorDataPair = getPassingPercolatorPeptides(
					parameters, pecanJobs, representativeJob, featureFiles, bigFeatureFile, bigPercolatorFile,
					bigPercolatorDecoyFile, bigPercolatorFiles, threshold, anyDDA);
			final Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides = percolatorDataPair.x;

			Logger.logLine("Identified "+passingPeptides.x.size()+" peptides across all files at a "+(threshold*100.0f)+"% FDR threshold.");

			boolean foundLibrary = checkLibraryForAvailability(libFile, percolatorDataPair.y);
			if (!foundLibrary) {
				Optional<PercolatorExecutionData> globalPercolatorFiles = Optional.ofNullable(featureFiles.size() == 1 ? null : bigPercolatorFiles);

				quantifySamples(progress, pecanJobs, libFile, outputFormat, alignBetweenFiles, parameters, percolatorDataPair.x, globalPercolatorFiles, integratePrecursors);
			}
			progress.update(percolatorDataPair.x.x.size()+" peptides identified at "+(threshold*100.0f)+"% FDR", 1.0f);
		} catch (IOException ioe) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ioe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error creating concatenated feature file");
			Logger.errorException(ie);
		}
	}

	private static Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Boolean> getPassingPercolatorPeptides(
			SearchParameters parameters, List<? extends SearchJobData> pecanJobs, SearchJobData representativeJob,
			ArrayList<File> featureFiles, File bigFeatureFile, File bigPercolatorFile, File bigPercolatorDecoyFile,
			PercolatorExecutionData bigPercolatorFiles, final float threshold, boolean anyDDA)
			throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides;
		boolean runningPercolator=true;
		if (featureFiles.size()==1) {
			Logger.logLine("Only one file, so no need to re-run Percolator.");
			// if there's only one file then don't need to re-run percolator
			passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(representativeJob.getPercolatorFiles().getPeptideOutputFile(), parameters, false);
			runningPercolator=false;
			
		} else if (parameters.isDoNotUseGlobalFDR()) {
			Logger.logLine("Warning, user asked to not use global FDR!");
			passingPeptides=getPeptidesWithoutGlobalFDR(pecanJobs, parameters).x;
			runningPercolator=false;
			
		} else if (bigPercolatorFile.exists()&&bigPercolatorFile.canRead()&&bigPercolatorDecoyFile.exists()&&bigPercolatorDecoyFile.canRead()) {
			Logger.logLine("Found previously run global Percolator.");
			// if we've already run percolator then don't need to re-run percolator
			passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(bigPercolatorFile, parameters, false);
			runningPercolator=false;
		} else {
			Logger.logLine("Running global Percolator analysis.");
			TableConcatenator.concatenatePINTables(featureFiles, bigFeatureFile, representativeJob.getPrimaryScoreName(), anyDDA);
			
			// delete if exists
			if (bigPercolatorFiles.getModelFile().exists()) {
				bigPercolatorFiles.getModelFile().delete();
			}
			int modelNumber = Integer.MAX_VALUE; // always use the last model (if reusing a model)
			

			if (parameters.isUsePercolator()) {
				Logger.logLine("Running JavaPot ("+(parameters.getPercolatorThreshold()*100f)+"%)");
				passingPeptides=JavaPotExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), bigPercolatorFiles, threshold, parameters.getAAConstants(), modelNumber);
			} else {
				Logger.logLine("Running mProphet ("+(parameters.getPercolatorTestThreshold()*100f)+"%)");
				MProphetExecutionData mprophetData=new MProphetExecutionData(bigPercolatorFiles);
				MProphetResult result=MProphetReiter.executeMProphetTSV(mprophetData, threshold, parameters.getAAConstants(), 1);
				passingPeptides=new Pair<ArrayList<PercolatorPeptide>, Float>(result.getPassingPeptides(), result.getPi0());
			}
		}
		
		Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Boolean> percolatorDataPair=new Pair<>(passingPeptides, runningPercolator);

		Logger.logLine("Identified "+percolatorDataPair.x.x.size()+" peptides across all files at a "+(threshold*100.0f)+"% FDR threshold.");
		return percolatorDataPair;
	}

	private static boolean checkLibraryForAvailability(File libFile, boolean runningPercolator) {
		boolean foundLibrary=false;
		if ((!runningPercolator)&&libFile.exists()&&libFile.canRead()) {
			// didn't have to run percolator, so check if we can read the lib file
			try {
				LibraryFile lib=new LibraryFile();
				lib.openFile(libFile);
				Logger.logLine("Found library file and tested for reading. It seems ok so proceeding with that file!");
				foundLibrary=true;
				
			} catch (Exception e) {
				Logger.logLine("Found library file and tested for reading. Reading failed, so overwriting!");
			}
		}
		return foundLibrary;
	}

	private static List<? extends SearchJobData> getProcessedJobs(List<? extends SearchJobData> origPecanJobs) {
		ArrayList<SearchJobData> processedJobs=new ArrayList<SearchJobData>();
		List<? extends SearchJobData> pecanJobsObj = Lists.newArrayList(origPecanJobs); // mutable copy
		// Sort files in alphabetical order for deterministic Percolator sampling
		Collections.sort(pecanJobsObj, (a, b) -> a.getOriginalDiaFileName().compareTo(b.getOriginalDiaFileName()));

		for (int i=0; i<pecanJobsObj.size(); i++) {
			SearchJobData job=pecanJobsObj.get(i);
			if (!job.hasBeenRun()) {
				Logger.logLine("Can't find a "+job.getSearchType()+" analysis of "+job.getOriginalDiaFileName()+", skipping extraction on that file.");
				continue;
			} else {
				processedJobs.add(job);
			}
		}

		return processedJobs;
	}

	private static void quantifySamples(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs,
										File libFile, OutputFormat outputFormat, boolean alignBetweenFiles, SearchParameters parameters,
										Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides,
										Optional<PercolatorExecutionData> globalPercolatorFiles, boolean integratePrecursors) {
		Optional<PeakLocationInferrerInterface> inferrer;
		if (alignBetweenFiles) {
			Logger.logLine("Inferring peak boundaries across files...");
			try {
				inferrer=Optional.of(EncyclopediaTwoPeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), pecanJobs, passingPeptides.x, parameters));
				Logger.logLine("...Finished peak inference.");
			} catch (Exception e) {
				Logger.errorLine("RT alignment between files failed! Perhaps this is to build a chromatogram library and not a quantitative experiment? Attempting to recover without alignment.");
				Logger.errorException(e);
				inferrer=Optional.empty();
			}
		} else {
			Logger.logLine("No RT alignment between files necessary.");
			inferrer=Optional.empty();
		}

		outputFormat.convert(progress, pecanJobs, libFile, passingPeptides, globalPercolatorFiles, inferrer, parameters, integratePrecursors);
	}

	private static Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Pair<ArrayList<PercolatorPeptide>, Float>> getPeptidesWithoutGlobalFDR(List<? extends SearchJobData> pecanJobs, SearchParameters parameters) {
		Pair<ArrayList<PercolatorPeptide>, Float> resultDecoyPeptides=new Pair<ArrayList<PercolatorPeptide>, Float>(new ArrayList<>(), -1f);
		HashMap<String, ScoredObject<PeptidePrecursor>> decoyMap=new HashMap<>();
		for (SearchJobData job : pecanJobs) {
			ArrayList<PercolatorPeptide> individualSamplePeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideDecoyFile(), parameters, false).x;
			for (PercolatorPeptide peptide : individualSamplePeptides) {
				ScoredObject<PeptidePrecursor> obj=decoyMap.get(peptide.getPeptideModSeq());
				if (obj==null||obj.x>peptide.getPosteriorErrorProb()) {
					decoyMap.put(peptide.getPeptideModSeq(), new ScoredObject<PeptidePrecursor>(peptide.getPosteriorErrorProb(), peptide));
				}
			}
		}
		for (ScoredObject<PeptidePrecursor> precursor : decoyMap.values()) {
			resultDecoyPeptides.x.add((PercolatorPeptide)precursor.y);
		}
		
		Pair<ArrayList<PercolatorPeptide>, Float> resultTargetPeptides=new Pair<ArrayList<PercolatorPeptide>, Float>(new ArrayList<>(), -1f);
		HashMap<String, ScoredObject<PeptidePrecursor>> targetMap=new HashMap<>();
		for (SearchJobData job : pecanJobs) {
			ArrayList<PercolatorPeptide> individualSamplePeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), parameters, false).x;
			for (PercolatorPeptide peptide : individualSamplePeptides) {
				ScoredObject<PeptidePrecursor> obj=targetMap.get(peptide.getPeptideModSeq());
				if (obj==null||obj.x>peptide.getPosteriorErrorProb()) {
					targetMap.put(peptide.getPeptideModSeq(), new ScoredObject<PeptidePrecursor>(peptide.getPosteriorErrorProb(), peptide));
				}
			}
		}
		for (ScoredObject<PeptidePrecursor> precursor : targetMap.values()) {
			resultTargetPeptides.x.add((PercolatorPeptide)precursor.y);
		}
		
		return new Pair<Pair<ArrayList<PercolatorPeptide>,Float>, Pair<ArrayList<PercolatorPeptide>,Float>>(resultTargetPeptides, resultDecoyPeptides);
	}
	
	static void convertBlib(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File blibFile, Optional<ArrayList<PercolatorPeptide>> passingPeptides, Optional<PeakLocationInferrerInterface> inferrer) {
		try {
			BlibFile blib=new BlibFile();
			blib.openFile();
			blib.setUserFile(blibFile);
			blib.dropIndices();
			int[] counterTotals=new int[] {0,0,0};
			
			File integrationFile=new File(blibFile.getAbsolutePath()+".integration.txt");
			PrintWriter integrationFileWriter=new PrintWriter(integrationFile, "UTF-8");
			integrationFileWriter.println("File Name\tPeptide Modified Sequence\tMin Start Time\tMax End Time\tPrecursor Charge\tPrecursorIsDecoy\tIon Count\tRetention Time Center\tTIC");

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<PercolatorPeptide> globalPassingPeptides;
				ArrayList<PercolatorPeptide> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), pecanJobs.get(i).getParameters(), false).x;
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get();
				} else {
					globalPassingPeptides=localPassingPeptides;
				}
				
				counterTotals=convertFileBlib(subProgress, job, globalPassingPeptides, localPassingPeptides, counterTotals, inferrer, integrationFileWriter, blib);
			}
			integrationFileWriter.flush();
			integrationFileWriter.close();

			blib.createIndices();
			blib.saveFile();
			blib.close();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	/**
	 * trims to quantifiable peptides! for loading into skyline!
	 */
	static int[] convertFileBlib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, int[] counterTotals, Optional<PeakLocationInferrerInterface> inferrer, PrintWriter integrationFileWriter, BlibFile blib) throws IOException, SQLException {
		final String diaFileName = job.getOriginalDiaFileName();

		Logger.logLine("Reading Percolator Results from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Reading Percolator Results", 0.0f);

		final StripeFileInterface stripeFile = job.getDiaFileReader();

		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		LibraryInterface library=null;
		if (job instanceof EncyclopediaJobData) {
			library=((EncyclopediaJobData)job).getLibrary();
		}
		//ArrayList<IntegratedLibraryEntry> libraryEntries=SearchFeatureReader.parseSearchFeatures(featureFile, globalPassingPeptides, localPassingPeptides, stripeFile, Optional.ofNullable((LibraryFile)null), job.getParameters());
		ArrayList<IntegratedLibraryEntry> libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, true, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
		stripeFile.close();
		
		for (IntegratedLibraryEntry entry : libraryEntries) {
			String peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(entry.getPeptideModSeq());
			integrationFileWriter.println(diaFileName +"\t"+peptideModSeq+"\t"+entry.getRtRange().getStart()/60f+"\t"+entry.getRtRange().getStop()/60f+"\t"+entry.getPrecursorCharge()+"\tFALSE\t"+entry.getIonCount()+"\t"+entry.getRetentionTime()/60f+"\t"+entry.getTIC());
		}
		integrationFileWriter.flush();
		
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : libraryEntries) {
			recasted.add(entry);
		}

		Logger.logLine("Writing Skyline BLIB from "+ diaFileName +"...");
		subProgress.update(diaFileName +": Writing Skyline BLIB", 0.99999f);

		counterTotals=blib.addLibrary(job, recasted, counterTotals[0], counterTotals[1], counterTotals[2]);
		subProgress.update(diaFileName +": Finished writing to Skyline BLIB at"+new Date().toString(), 1.0f);
		return counterTotals;
	}

	static void convertElib(ProgressIndicator progress, SearchJobData pecanJob, File elibFile, SearchParameters parameters, boolean integratePrecursors) {
		ArrayList<SearchJobData> jobs=new ArrayList<>();
		jobs.add(pecanJob);

		convertElib(progress, jobs, elibFile, Optional.empty(), Optional.empty(), Optional.empty(), parameters, integratePrecursors);
	}

	static void convertElib(ProgressIndicator progress, List<? extends SearchJobData> pecanJobs, File elibFile, Optional<Pair<ArrayList<PercolatorPeptide>, Float>> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, Optional<PeakLocationInferrerInterface> inferrer, SearchParameters parameters, boolean integratePrecursors) {
		try {
			long convertStart=System.currentTimeMillis();
			logScribeTwoPerf("result_library_start jobs="+pecanJobs.size()+" output="+elibFile.getName()+" integratePrecursors="+integratePrecursors);
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();

			float increment=1.0f/pecanJobs.size();
			for (int i=0; i<pecanJobs.size(); i++) {
				SearchJobData job=pecanJobs.get(i);
				if (!job.hasBeenRun()) {
					continue;
				}
				ProgressIndicator subProgress=new SubProgressIndicator(progress, increment);
				
				ArrayList<PercolatorPeptide> globalPassingPeptides;
				long localPeptideReadStart=System.currentTimeMillis();
				Pair<ArrayList<PercolatorPeptide>, Float> localPassingPeptides=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), pecanJobs.get(i).getParameters(), false);
				logScribeTwoPerf("result_library_local_peptides file="+job.getOriginalDiaFileName()+" peptides="+localPassingPeptides.x.size()
						+" elapsedMs="+(System.currentTimeMillis()-localPeptideReadStart));
				if (passingPeptides.isPresent()) {
					globalPassingPeptides=passingPeptides.get().x;
				} else {
					globalPassingPeptides=localPassingPeptides.x;
				}

				Logger.logLine(job.getOriginalDiaFileName()+": Number of global peptides: "+globalPassingPeptides.size()+" vs local peptides: "+localPassingPeptides.x.size());
				
				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides.x, inferrer, elib, pecanJobs.size()>1, integratePrecursors);

				if ((!globalPercolatorFiles.isPresent())) {
					if (job.hasBeenRun()) {
						long targetDecoyStart=System.currentTimeMillis();
						Pair<ArrayList<PercolatorPeptide>, Float> targets=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), parameters, true);
						Pair<ArrayList<PercolatorPeptide>, Float> decoys=PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideDecoyFile(), parameters, true);
						Logger.logLine("Writing local target/decoy peptides: "+targets.x.size()+"/"+decoys.x.size()+", pi0: "+targets.y);
						elib.addTargetDecoyPeptides(targets.x, decoys.x);
						elib.addMetadata("pi0", Float.toString(targets.y));
						elib.addProteinsFromPercolator(targets.x);
						elib.addProteinsFromPercolator(decoys.x);
						logScribeTwoPerf("result_library_target_decoy_peptides file="+job.getOriginalDiaFileName()+" targets="+targets.x.size()
								+" decoys="+decoys.x.size()+" elapsedMs="+(System.currentTimeMillis()-targetDecoyStart));
						
						long proteinInferenceStart=System.currentTimeMillis();
						Pair<ArrayList<PercolatorPeptide>, Float> proteinInferenceTargets=PercolatorReader.getScoredPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), parameters.getAAConstants(), true);
						Pair<ArrayList<PercolatorPeptide>, Float> proteinInferenceDecoys=PercolatorReader.getScoredPeptidesFromTSV(job.getPercolatorFiles().getPeptideDecoyFile(), parameters.getAAConstants(), true);
						PickedProteinFDRResult proteinFDR=ParsimonyProteinGrouper.calculatePickedProteinFDR(proteinInferenceTargets.x, proteinInferenceDecoys.x,
								parameters.getPercolatorProteinThreshold(), parameters.getAAConstants(), parameters.getProteinGroupScoreMode());
						logProteinInferenceSummary("local", proteinFDR);
						Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=
								new Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>>(proteinFDR.getAcceptedTargets(), proteinFDR.getAcceptedDecoys());
						targetDecoyProteins=new Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>>(
								filterProteinGroupsWithPassingPeptideEvidence(targetDecoyProteins.x, targets.x),
								filterProteinGroupsWithPassingPeptideEvidence(targetDecoyProteins.y, decoys.x));
						Logger.logLine("Writing local target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
						logEntrapmentSummaries("local", targets.x, decoys.x, targetDecoyProteins.x, targetDecoyProteins.y);
						elib.addTargetDecoyProteins(job.getOriginalDiaFileName(), targetDecoyProteins.x, targetDecoyProteins.y);
						logScribeTwoPerf("result_library_protein_inference file="+job.getOriginalDiaFileName()+" targets="+targetDecoyProteins.x.size()
								+" decoys="+targetDecoyProteins.y.size()+" elapsedMs="+(System.currentTimeMillis()-proteinInferenceStart));

						job.getPercolatorFiles()
								.getPercolatorExecutableVersion()
								.ifPresent((ThrowingConsumer<String>) version -> {
									elib.addMetadata(LibraryFile.PERCOLATOR_VERSION, version);
								});
					}
				}
				
				subProgress.update("Wrote "+globalPassingPeptides.size()+" peptides identified at "+(job.getParameters().getPercolatorThreshold()*100.0f)+"% FDR", 1.0f);
			}
			

			ArrayList<PercolatorProteinGroup> proteins=null;
			if (globalPercolatorFiles.isPresent()) {
				long globalPercolatorStart=System.currentTimeMillis();
				proteins = writePercolatorToElib(elib, globalPercolatorFiles.get(), pecanJobs, parameters);
				logScribeTwoPerf("result_library_global_percolator proteins="+(proteins==null?0:proteins.size())+" elapsedMs="
						+(System.currentTimeMillis()-globalPercolatorStart));
			}

			long metadataStart=System.currentTimeMillis();
			writeElibMetadata(elib, pecanJobs, parameters, inferrer.isPresent());
			logScribeTwoPerf("result_library_metadata elapsedMs="+(System.currentTimeMillis()-metadataStart));

			long createIndicesStart=System.currentTimeMillis();
			elib.createIndices();
			logScribeTwoPerf("result_library_create_indices elapsedMs="+(System.currentTimeMillis()-createIndicesStart));
			long saveStart=System.currentTimeMillis();
			elib.saveAsFile(elibFile);
			logScribeTwoPerf("result_library_save elapsedMs="+(System.currentTimeMillis()-saveStart));

			if (proteins!=null) {
				if (inferrer.isPresent()) {
					try {
						long reportStart=System.currentTimeMillis();
						LibraryReportExtractor.extractMatrix(elib, proteins, parameters.isNormalizeByTIC());
						logScribeTwoPerf("result_library_extract_matrix elapsedMs="+(System.currentTimeMillis()-reportStart));
					} catch (DataFormatException e) {
						Logger.errorException(e);
					}
				} else {
					Logger.errorLine("Only exporting report for a single search, so skipping building quantitative tables.");
				}
			}

			elib.close();
			logScribeTwoPerf("result_library_complete elapsedMs="+(System.currentTimeMillis()-convertStart));
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	/**
	 * Does not limit to quantifiable! Reports all potential peaks!
	 */
	static void convertFileElib(ProgressIndicator subProgress, SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPeptides, ArrayList<PercolatorPeptide> localPassingPeptides, Optional<PeakLocationInferrerInterface> inferrer, LibraryFile elib, boolean combineJobs, boolean integratePrecursors) throws IOException, SQLException {
		long fileStart=System.currentTimeMillis();
		String diaFileName=job.getOriginalDiaFileName();
		Logger.logLine("Reading Percolator Results from "+diaFileName+"...");
		subProgress.update(diaFileName+": Reading Percolator Results", 0.0f);

		long readerStart=System.currentTimeMillis();
		final StripeFileInterface stripeFile = job.getDiaFileReader();
		logScribeTwoPerf("result_library_reader_open file="+diaFileName+" elapsedMs="+(System.currentTimeMillis()-readerStart));

		Logger.logLine("Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides from "+diaFileName+"...");
		subProgress.update(diaFileName+": Extracting Spectral Data for "+localPassingPeptides.size()+" Peptides", 0.00001f);

		long ticStart=System.currentTimeMillis();
		elib.addTIC(stripeFile);
		logScribeTwoPerf("result_library_tic file="+diaFileName+" elapsedMs="+(System.currentTimeMillis()-ticStart));

		inferrer.ifPresent(inf -> elib.addRtAlignment(job, inf));
		boolean hasLocalizationData=!combineJobs&&(job instanceof ThesaurusJobData||job instanceof VariantXCorDIAJobData);
		boolean useStreamingDDAEntries=integratePrecursors&&!hasLocalizationData;

		ArrayList<IntegratedLibraryEntry> libraryEntries=null;
		int libraryEntryCount;
		boolean integratedEntriesWritten=false;
		if (job instanceof QuantitativeSearchJobData&&(!integratePrecursors)) {
			LibraryInterface library=null;
			if (job instanceof EncyclopediaJobData) {
				library=((EncyclopediaJobData)job).getLibrary();
			}
			
			long quantStart=System.currentTimeMillis();
			libraryEntries=PeptideQuantExtractor.parseSearchFeatures(subProgress, job, false, globalPassingPeptides, localPassingPeptides, inferrer, stripeFile, library, job.getParameters());
			libraryEntryCount=libraryEntries.size();
			logScribeTwoPerf("result_library_fragment_quant file="+diaFileName+" entries="+libraryEntries.size()+" elapsedMs="
					+(System.currentTimeMillis()-quantStart));
		} else {
			long targetPsmStart=System.currentTimeMillis();
			HashMap<String, PSMData> targetPSMs=PeptideQuantExtractor.findTargetPSMData(job, globalPassingPeptides, localPassingPeptides, inferrer, job.getParameters());
			logScribeTwoPerf("result_library_target_psms file="+diaFileName+" psms="+targetPSMs.size()+" elapsedMs="
					+(System.currentTimeMillis()-targetPsmStart));
			long ddaQuantStart=System.currentTimeMillis();
			if (useStreamingDDAEntries) {
				try (final LibraryFile.DDAIntegratedEntryWriter writer=elib.openDDAIntegratedEntryWriter(!integratePrecursors, inferrer, job.getParameters())) {
					libraryEntryCount=DDAPrecursorIntegrator.integrateSearch(subProgress, new ArrayList<PSMData>(targetPSMs.values()), stripeFile, job.getParameters(),
							new DDAPrecursorIntegrator.IntegratedEntryConsumer() {
								@Override
								public void accept(ArrayList<IntegratedLibraryEntry> entries) throws IOException, SQLException {
									writer.addEntries(entries);
								}
							});
					writer.finish();
				}
				integratedEntriesWritten=true;
			} else {
				libraryEntries=DDAPrecursorIntegrator.integrateSearch(subProgress, new ArrayList<PSMData>(targetPSMs.values()), stripeFile, job.getParameters());
				libraryEntryCount=libraryEntries.size();
			}
			logScribeTwoPerf("result_library_dda_quant file="+diaFileName+" entries="+libraryEntryCount+" elapsedMs="
					+(System.currentTimeMillis()-ddaQuantStart));
		}
		stripeFile.close();
		
		Logger.logLine("Writing Encyclopedia ELIB from "+diaFileName+" ("+libraryEntryCount+" entries)...");
		subProgress.update(diaFileName+": Writing Encyclopedia ELIB", 0.99999f);
		
		Optional<HashMap<String, ModificationLocalizationData>> localizationData;
		if (hasLocalizationData&&job instanceof ThesaurusJobData) {
			Logger.logLine("Reading localization data from disk...");
			localizationData=Optional.of(LocalizationDataToTSVConsumer.readLocalizationFile(((ThesaurusJobData)job).getLocalizationFile(), globalPassingPeptides, job.getParameters()));
		} else if (hasLocalizationData&&job instanceof VariantXCorDIAJobData) {
			Logger.logLine("Reading localization data from disk...");
			localizationData=Optional.of(LocalizationDataToTSVConsumer.readLocalizationFile(((VariantXCorDIAJobData)job).getLocalizationFile(), globalPassingPeptides, job.getParameters()));
		} else {
			localizationData=Optional.empty();
		}

		long entryWriteStart=System.currentTimeMillis();
		if (!integratedEntriesWritten) {
			elib.addIntegratedEntries(!integratePrecursors, libraryEntries, inferrer, localizationData, job.getParameters());
		}
		logScribeTwoPerf("result_library_entries_written file="+diaFileName+" entries="+libraryEntryCount+" streamed="+integratedEntriesWritten+" elapsedMs="
				+(System.currentTimeMillis()-entryWriteStart));
		

		Logger.logLine("Finished writing to Encyclopedia ELIB at "+new Date().toString());
		subProgress.update(diaFileName+": Finished writing to Encyclopedia ELIB at "+new Date().toString(), 1.0f);
		logScribeTwoPerf("result_library_file_complete file="+diaFileName+" elapsedMs="+(System.currentTimeMillis()-fileStart));
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata. If the
	 *                                global results file doesn't exist this will be ignored and the calculation will
	 *                                fall back to use {@code jobs}.
	 * @param jobs Ignored, unless global Percolator results don't exist, in which case the passing peptides are read
	 *             directly from these jobs, without global FDR control.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, List<? extends SearchJobData> jobs, SearchParameters parameters) throws IOException, SQLException {
		return writePercolatorToElib(elib, percolatorExecutionData, Optional.of(jobs), parameters);
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, SearchParameters parameters) throws IOException, SQLException {
		return writePercolatorToElib(elib, percolatorExecutionData, Optional.empty(), parameters);
	}

	/**
	 * Read the set of passing peptides from Percolator and write them to the ELIB with associated metadata.
	 * Perform protein inference and write protein scores/q-values/PEPs to the ELIB.
	 *
	 * @param elib The (open) ELIB where results will be written.
	 * @param percolatorExecutionData Used to read the list of passing peptides and associated scores/metadata. If the
	 *                                global results file doesn't exist this will be ignored and the calculation will
	 *                                fall back to use {@code jobs}.
	 * @param jobs Ignored, unless global Percolator results don't exist, in which case the passing peptides are read
	 *             directly from these jobs, without global FDR control. An exception will be raised if this fallback
	 *             is necessary but {@code jobs} is not present.
	 * @return The inferred set of protein groups.
	 */
	private static ArrayList<PercolatorProteinGroup> writePercolatorToElib(LibraryFile elib, PercolatorExecutionData percolatorExecutionData, Optional<List<? extends SearchJobData>> jobs, SearchParameters parameters) throws IOException, SQLException {
		Pair<ArrayList<PercolatorPeptide>, Float> targets=null;
		Pair<ArrayList<PercolatorPeptide>, Float> decoys=null;
		Pair<ArrayList<PercolatorPeptide>, Float> proteinInferenceTargets=null;
		Pair<ArrayList<PercolatorPeptide>, Float> proteinInferenceDecoys=null;
		if (percolatorExecutionData.getPeptideOutputFile().exists()) {
			targets=PercolatorReader.getPassingPeptidesFromTSV(percolatorExecutionData.getPeptideOutputFile(), parameters, true);
			decoys=PercolatorReader.getPassingPeptidesFromTSV(percolatorExecutionData.getPeptideDecoyFile(), parameters, true);
			proteinInferenceTargets=PercolatorReader.getScoredPeptidesFromTSV(percolatorExecutionData.getPeptideOutputFile(), parameters.getAAConstants(), true);
			proteinInferenceDecoys=PercolatorReader.getScoredPeptidesFromTSV(percolatorExecutionData.getPeptideDecoyFile(), parameters.getAAConstants(), true);
		} else if (jobs.isPresent()) {
			Pair<Pair<ArrayList<PercolatorPeptide>, Float>, Pair<ArrayList<PercolatorPeptide>, Float>> withoutFDR=getPeptidesWithoutGlobalFDR(jobs.get(), parameters);
			targets=withoutFDR.x;
			decoys=withoutFDR.y;
			proteinInferenceTargets=targets;
			proteinInferenceDecoys=decoys;
			Logger.logLine("Broad protein inference TSV evidence unavailable; using no-global-FDR thresholded peptide lists for protein inference.");
		} else {
			throw new IllegalStateException("Unable to get passing peptides: no global Percolator results file or individual jobs!");
		}

		Logger.logLine("Writing global target/decoy peptides: "+targets.x.size()+"/"+decoys.x.size()+", pi0: "+targets.y);
		elib.addTargetDecoyPeptides(targets.x, decoys.x);
		elib.addMetadata("pi0", Float.toString(targets.y));
		elib.addProteinsFromPercolator(targets.x);
		elib.addProteinsFromPercolator(decoys.x);

		PickedProteinFDRResult proteinFDR=ParsimonyProteinGrouper.calculatePickedProteinFDR(proteinInferenceTargets.x, proteinInferenceDecoys.x,
				parameters.getPercolatorProteinThreshold(), parameters.getAAConstants(), parameters.getProteinGroupScoreMode());
		logProteinInferenceSummary("global", proteinFDR);
		Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=
				new Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>>(proteinFDR.getAcceptedTargets(), proteinFDR.getAcceptedDecoys());
		targetDecoyProteins=new Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>>(
				filterProteinGroupsWithPassingPeptideEvidence(targetDecoyProteins.x, targets.x),
				filterProteinGroupsWithPassingPeptideEvidence(targetDecoyProteins.y, decoys.x));

		Logger.logLine("Writing global target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
		logEntrapmentSummaries("global", targets.x, decoys.x, targetDecoyProteins.x, targetDecoyProteins.y);
		elib.addTargetDecoyProteins(GLOBAL_NAME, targetDecoyProteins.x, targetDecoyProteins.y);

		percolatorExecutionData
				.getPercolatorExecutableVersion()
				.ifPresent((ThrowingConsumer<String>) version -> {
					elib.addMetadata(LibraryFile.PERCOLATOR_VERSION, version);
				});

		return targetDecoyProteins.x;
	}

	private static void logProteinInferenceSummary(String scope, PickedProteinFDRResult result) {
		Logger.logLine("Protein inference "+scope+" evidence: target peptides="+result.getTargetPeptides()+", decoy peptides="+result.getDecoyPeptides()
				+", target groups="+result.getGroupedTargetGroups()+", decoy groups="+result.getGroupedDecoyGroups()
				+", paired keys="+result.getPairedKeys()+", target-only keys="+result.getTargetOnlyKeys()+", decoy-only keys="+result.getDecoyOnlyKeys()
				+", accepted target groups="+result.getAcceptedTargetGroups()+", accepted decoy groups="+result.getAcceptedDecoyGroups());
	}

	static ArrayList<PercolatorProteinGroup> filterProteinGroupsWithPassingPeptideEvidence(ArrayList<PercolatorProteinGroup> proteinGroups, ArrayList<PercolatorPeptide> passingPeptides) {
		HashSet<String> accessionsWithPassingPeptides=new HashSet<String>();
		for (PercolatorPeptide peptide : passingPeptides) {
			accessionsWithPassingPeptides.addAll(peptide.getAccessions());
		}

		ArrayList<PercolatorProteinGroup> filteredProteinGroups=new ArrayList<PercolatorProteinGroup>();
		for (PercolatorProteinGroup proteinGroup : proteinGroups) {
			for (String accession : proteinGroup.getEquivalentAccessions()) {
				if (accessionsWithPassingPeptides.contains(accession)) {
					filteredProteinGroups.add(proteinGroup);
					break;
				}
			}
		}
		return filteredProteinGroups;
	}

	public static void logEntrapmentSummaries(String scope, ArrayList<PercolatorPeptide> targetPeptides, ArrayList<PercolatorPeptide> decoyPeptides,
			ArrayList<? extends ProteinGroupInterface> targetProteins, ArrayList<? extends ProteinGroupInterface> decoyProteins) {
		EntrapmentSummary peptideSummary=getPeptideEntrapmentSummary(targetPeptides, decoyPeptides);
		EntrapmentSummary proteinSummary=getProteinEntrapmentSummary(targetProteins, decoyProteins);
		if (!peptideSummary.hasEntrapment()&&!proteinSummary.hasEntrapment()) return;

		Logger.logLine("Entrapment "+scope+" peptides: targets="+peptideSummary.targets+", decoys="+peptideSummary.decoys
				+", target-entrapments="+peptideSummary.targetEntrapments+", decoy-entrapments="+peptideSummary.decoyEntrapments
				+", lower-bound FDR="+formatEntrapmentFDR(peptideSummary.getLowerBoundFDR())
				+", upper-bound FDR="+formatEntrapmentFDR(peptideSummary.getUpperBoundFDR())
				+", assessment="+peptideSummary.getAssessment()
				);
		Logger.logLine("Entrapment "+scope+" proteins: targets="+proteinSummary.targets+", decoys="+proteinSummary.decoys
				+", target-entrapments="+proteinSummary.targetEntrapments+", decoy-entrapments="+proteinSummary.decoyEntrapments
				+", lower-bound FDR="+formatEntrapmentFDR(proteinSummary.getLowerBoundFDR())
				+", upper-bound FDR="+formatEntrapmentFDR(proteinSummary.getUpperBoundFDR())
				+", assessment="+proteinSummary.getAssessment()
				);
	}

	static EntrapmentSummary getPeptideEntrapmentSummary(ArrayList<PercolatorPeptide> targetPeptides, ArrayList<PercolatorPeptide> decoyPeptides) {
		EntrapmentSummary summary=new EntrapmentSummary(targetPeptides.size(), decoyPeptides.size());
		for (PercolatorPeptide peptide : targetPeptides) {
			if (hasEntrapmentAccession(peptide.getAccessions())) summary.targetEntrapments++;
		}
		for (PercolatorPeptide peptide : decoyPeptides) {
			if (hasEntrapmentAccession(peptide.getAccessions())) summary.decoyEntrapments++;
		}
		return summary;
	}

	static EntrapmentSummary getProteinEntrapmentSummary(ArrayList<? extends ProteinGroupInterface> targetProteins, ArrayList<? extends ProteinGroupInterface> decoyProteins) {
		EntrapmentSummary summary=new EntrapmentSummary(targetProteins.size(), decoyProteins.size());
		for (ProteinGroupInterface protein : targetProteins) {
			if (hasEntrapmentAccession(protein.getEquivalentAccessions())) summary.targetEntrapments++;
		}
		for (ProteinGroupInterface protein : decoyProteins) {
			if (hasEntrapmentAccession(protein.getEquivalentAccessions())) summary.decoyEntrapments++;
		}
		return summary;
	}

	static boolean hasEntrapmentAccession(Iterable<String> accessions) {
		for (String accession : accessions) {
			if (accession.endsWith("_p_target")) return true;
		}
		return false;
	}

	static String formatEntrapmentFDR(float fdr) {
		return String.format(Locale.US, "%.2f%%", 100.0f*fdr);
	}

	static class EntrapmentSummary {
		final int targets;
		final int decoys;
		int targetEntrapments;
		int decoyEntrapments;

		EntrapmentSummary(int targets, int decoys) {
			this.targets=targets;
			this.decoys=decoys;
		}

		boolean hasEntrapment() {
			return targetEntrapments>0||decoyEntrapments>0;
		}
		
		String getAssessment() {
			if (getLowerBoundFDR()>0.01f) {
				return "anti-conservative";
			} else if (getUpperBoundFDR()<0.01f) {
				return "overly-conservative";
			} else {
				return "well-controlled";
			}
		}

		float getLowerBoundFDR() {
			return targetEntrapments/(float)Math.max(targets, 1);
		}
		
		float getUpperBoundFDR() {
			return (2.0f*targetEntrapments)/(float)Math.max(targets+targetEntrapments, 1);
		}
	}

	private static void writeElibMetadata(LibraryFile elib, List<? extends SearchJobData> jobs, SearchParameters parameters, boolean align) throws IOException, SQLException {
		final HashMap<String, String> parameterMap = parameters.toParameterMap();
		parameterMap.put("RT align between samples", Boolean.toString(align));
		for (int i = 0; i < jobs.size(); i++) {
			final SearchJobData job = jobs.get(i);
			parameterMap.put(job.getOriginalDiaFileName() + " search type", job.getSearchType());
			if (job instanceof EncyclopediaJobData) {
				parameterMap.put(job.getOriginalDiaFileName() + " library", ((EncyclopediaJobData) job).getLibrary().getName());
			} else if (job instanceof PecanJobData) {
				parameterMap.put(job.getOriginalDiaFileName() + " fasta", ((PecanJobData) job).getFastaFile().getName());
				parameterMap.put(job.getOriginalDiaFileName() + " used narrow target list", Boolean.toString(((PecanJobData) job).getTargetList().isPresent()));
			} else if (job instanceof XCorDIAJobData) {
				Optional<LibraryInterface> maybeLibrary = ((XCorDIAJobData) job).getLibrary();
				if (maybeLibrary.isPresent()) {
					parameterMap.put(job.getOriginalDiaFileName() + " library", maybeLibrary.get().getName());
				}
				parameterMap.put(job.getOriginalDiaFileName() + " fasta", ((XCorDIAJobData) job).getFastaFile().getName());
				parameterMap.put(job.getOriginalDiaFileName() + " used narrow target list", Boolean.toString(((XCorDIAJobData) job).getTargetList().isPresent()));
			}
		}
		elib.addMetadata(parameterMap);

		elib.setSources(jobs);
	}

	private static void convertAlib(ProgressIndicator progress, List<? extends SearchJobData> jobs, File outputFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, Optional<PercolatorExecutionData> globalPercolatorFiles, PeakLocationInferrerInterface inferrer, SearchParameters parameters) {
		if (Objects.requireNonNull(jobs, "No jobs provided").isEmpty()) {
			throw new IllegalArgumentException("No jobs provided");
		}

		if (!parameters.isQuantifySameFragmentsAcrossSamples()) {
			throw new IllegalArgumentException("Unable to export alignment-only library without -quantifyAcrossSamples!");
		}

		if (!globalPercolatorFiles.isPresent()) {
			if (jobs.size() == 1) {
				globalPercolatorFiles = Optional.of(jobs.iterator().next().getPercolatorFiles());
			} else {
				throw new IllegalArgumentException("Global percolator files must be provided for more than one job!");
			}
		}

		try {
			final LibraryFile elib = new LibraryFile();
			try {
				elib.openFile();
				elib.dropIndices();

				final PercolatorExecutionData percolatorExecutionData = globalPercolatorFiles.get();
				if (!percolatorExecutionData.getPeptideOutputFile().exists()) {
					throw new IllegalArgumentException("Could not read Percolator results!", new FileNotFoundException(percolatorExecutionData.getPeptideOutputFile().getAbsolutePath()));
				}

				// Perform protein inference and write peptide/protein scores, metadata to ELIB
				writePercolatorToElib(elib, percolatorExecutionData, parameters);

				// Now compute and write the set of entries that to capture the alignment/transition refinement
				elib.addEntries(getAlignmentEntries(passingPeptides.x, jobs, inferrer, parameters));

				// Write information about each job to the ELIB
				float increment = 1.0f / jobs.size();
				for (int i = 0; i < jobs.size(); i++) {
					final SearchJobData job = jobs.get(i);

					final ProgressIndicator subProgress = new SubProgressIndicator(progress, increment);

					if (!job.hasBeenRun()) {
						final String msg = "Skipping incomplete job: " + job.getOriginalDiaFileName();
						subProgress.update(msg, 1f);
						Logger.errorLine(msg);
						continue;
					}

					elib.addRtAlignment(job, getRawAlignmentPoints(job, inferrer));

					elib.addTIC(job.getDiaFileReader());

					subProgress.update("Done writing job " + job.getOriginalDiaFileName(), 1f);
				}

				writeElibMetadata(elib, jobs, parameters, true); // align is required for ALIB

				elib.createIndices();
				elib.saveAsFile(outputFile);
			} finally {
				elib.close();
			}
		} catch (IOException | SQLException ioe) {
			Logger.errorLine("Error creating ELIB file");
			Logger.errorException(ioe);
			throw new EncyclopediaException("Error creating ELIB file", ioe);
		}
	}

	/**
	 * Create a special set of alignment points meant for writing to an ALIB, in order to capture the precise
	 * RT alignment mapping in each file. This is only possible in specific circumstances, but these should be
	 * satisfied for any run producing an ALIB (the inferrer, alignment, and warper all must be of the expected
	 * types, or this method falls back to match the normal set of points written to an ELIB).
	 *
	 * @see #readInferrer(LibraryFile, List, SearchParameters)
	 */
	private static List<AlignmentDataPoint> getRawAlignmentPoints(SearchJobData job, PeakLocationInferrerInterface inferrer) {
		if (inferrer instanceof SimplePeakLocationInferrer) {
			final SimplePeakLocationInferrer spli = (SimplePeakLocationInferrer) inferrer;
			final RetentionTimeAlignmentInterface alignment = spli.getAlignment(job);

			if (alignment instanceof AbstractRetentionTimeFilter) {
				AbstractRetentionTimeFilter artf = (AbstractRetentionTimeFilter) alignment;
				final edu.washington.gs.maccoss.encyclopedia.utils.math.Function rtWarper = artf.getRtWarper();

				if (rtWarper instanceof LinearInterpolatedFunction) {
					// Concatenate the peptides' points with the full set of knots from the function.
					// This should provide all required information while keeping all segment slopes the same.
					return Stream.concat(
							inferrer.getAlignmentData(job).stream(),
							rtWarper.getKnots().stream()
									.map(xy -> AlignmentDataPoint.of(
											(float) xy.x,
											Float.NEGATIVE_INFINITY, // work around NOT NULL constraint
											(float) xy.y,
											Float.NEGATIVE_INFINITY, // work around NOT NULL constraint
											Float.NEGATIVE_INFINITY, // work around NOT NULL constraint
											false,
											null
									))
					).collect(Collectors.toList());
				}
			}
		}

		// fallback -- not ideal because it only uses the input points (along with their predictions though)

		Logger.errorLine("Unable to get exact RT alignment mapping! Falling back to peptide retentiontimes only.");
		return inferrer.getAlignmentData(job);
	}

	/**
	 * Compute the set of entries to be written to an alignment-only ELIB (ALIB).
	 *
	 * This consists of one entry per item of {@code passingPeptides}, with the corresponding aligned RT.
	 * Each entry's set of peaks will be determined by the quantitative peaks returned from {@code inferrer}
	 * for that peptide.
	 *
	 * @see OutputFormat#ALIB
	 *
	 * @param passingPeptides The set of passing peptide IDs from Percolator
	 * @param jobs
	 * @param inferrer The retention time alignment and transition refinement results for the experiment. Importantly,
	 *                 all entries will include only the quantitative ions from this inferrer.
	 *
	 * @return A set of entries suitable for insertion in the ELIB (ALIB) file.
	 */
	private static ArrayList<LibraryEntry> getAlignmentEntries(List<? extends PercolatorPeptide> passingPeptides, List<? extends SearchJobData> jobs, PeakLocationInferrerInterface inferrer, SearchParameters parameters) {
		return passingPeptides.stream()
				.map(p -> toAlignmentEntry(p, jobs, inferrer, parameters))
				.collect(Collectors.toCollection(ArrayList::new));
	}

	private static LibraryEntry toAlignmentEntry(PercolatorPeptide peptide, List<? extends SearchJobData> jobs, PeakLocationInferrerInterface inferrer, SearchParameters parameters) {
		float warpedRTInSec;
		try {
			// We want the aligned ("seed") RT, not the time in any specific sample, so we pass a "bogus" job. TODO: BIG RISK (NPE, interface abuse)
			warpedRTInSec = inferrer.getWarpedRTInSec(null, peptide.getPeptideModSeq());
		} catch (NullPointerException e) {
			warpedRTInSec = -1f;
		}

		final double[] quantIons = inferrer.getTopNBestIons(peptide.getPeptideModSeq(), peptide.getPrecursorCharge());

		// Look up the original entry from the search library.
		// Be sure we get it from the correct job's library, as they aren't guaranteed to be the same!
		final String file = peptide.getFile();

		LibraryEntry entry = null;

		final Set<String> checkedLibraries = Sets.newHashSet();
		for (SearchJobData job : jobs) {
			if (!Objects.equals(file, job.getOriginalDiaFileName())) {
				// Skip jobs that don't match this peptide ID
				continue;
			}

			if (!(job instanceof EncyclopediaJobData)) {
				// Skip jobs without a library
				continue;
			}

			final LibraryInterface library = ((EncyclopediaJobData) job).getLibrary();
			if (!checkedLibraries.add(library.getName())) {
				// Already checked this library; skip it
				continue;
			}

			final ArrayList<LibraryEntry> entries;
			try {
				entries = library.getEntries(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), false);
			} catch (IOException | SQLException | DataFormatException e) {
				// Unable to read entries from library! We don't want to excessively log, so just move on.
				continue;
			}

			// Take the highest-scoring peptide (assumes higher scores are better)
			final Optional<LibraryEntry> bestEntry = entries.stream().max(Comparator.comparing(LibraryEntry::getScore));

			if (!bestEntry.isPresent()) {
				// Try a different job's library
				continue;
			}

			entry = bestEntry.get();
			break;
		}

		if (null == entry) {
			throw new IllegalStateException("Unable to find entry in original search library for " + peptide.getPsmID());
		}

		final double[] entryMasses = entry.getMassArray();
		final float[] entryIntensities = entry.getIntensityArray();
		final float[] entryCorrelations = entry.getCorrelationArray();

		final MassTolerance tol = parameters.getLibraryFragmentTolerance();

		// will initialize to all false
		final boolean[] entryQuantifiedIons = new boolean[entryMasses.length];

		// Final arrays may differ from the entry (if we have to insert additional quant ions)
		final double[] masses;
		final float[] intensities, correlations;
		final boolean[] quantifiedIons;

		// Quant ions array is sometimes null; treat as though it's empty
		if (null != quantIons) {
			// Could be faster if we avoid quadratic search but we can't be certain
			// that either array is sorted, and the # of quant ions should be small.

			final TDoubleSet quantIonSet = new TDoubleHashSet(quantIons);

			for (int i = 0; i < entryMasses.length; i++) {
				for(TDoubleIterator iterator = quantIonSet.iterator(); iterator.hasNext();) {
					final double quantIon = iterator.next();

					if (tol.equals(entryMasses[i], quantIon)) {
						entryQuantifiedIons[i] = true;
						iterator.remove(); // don't use this quant ion again, we already found it

						break;
					}
				}
			}

			if (quantIonSet.isEmpty()) {
				masses = entryMasses;
				intensities = entryIntensities;
				correlations = entryCorrelations;
				quantifiedIons = entryQuantifiedIons;
			} else {
				// Any quant ions not found in the entry must be inserted.

				final int len = entryMasses.length + quantIonSet.size();
				masses = new double[len];
				intensities = new float[len];
				correlations = new float[len];
				quantifiedIons = new boolean[len];

				final double[] toInsert = quantIonSet.toArray();
				Arrays.sort(toInsert);

				// Just assume masses are sorted -- will do weird stuff if not but _should_ work.
				// i -- index in entry
				// j -- index in toInsert
				// k -- index in final arrays
				int j = 0, k = 0;
				for (int i = 0; i < entryMasses.length; i++) {
					for (; j < toInsert.length && toInsert[j] < entryMasses[i]; j++) {
						// Time to insert a mass
						masses[k] = toInsert[j];
						intensities[k] = 0f;
						correlations[k] = 0f;
						quantifiedIons[k] = true;

						k += 1; // move to next insertion location
					}

					// Now the next mass to insert is greater than the current entry mass
					masses[k] = entryMasses[i];
					intensities[k] = entryIntensities[i];
					correlations[k] = entryCorrelations[i];
					quantifiedIons[k] = entryQuantifiedIons[i];

					k += 1; // move to next insertion location
				}

				// Insert remaining masses at end
				for (; j < toInsert.length; j++) {
					masses[k] = toInsert[j];
					intensities[k] = 0f;
					correlations[k] = 0f;
					quantifiedIons[k] = true;

					k += 1; // move to next insertion location
				}
			}

			// Sanity checks for quant ions
			if (quantIons.length != General.sum(quantifiedIons)) {
				throw new IllegalStateException(
						"Unable to locate all quantitative ions for "
								+ peptide.getPeptideModSeq()
								+ ": expected "
								+ quantIons.length
								+ " found "
								+ General.sum(quantifiedIons)
				);
			}
		} else {
			masses = entryMasses;
			intensities = entryIntensities;
			correlations = entryCorrelations;
			quantifiedIons = entryQuantifiedIons;
		}

		return new LibraryEntry(
				GLOBAL_NAME,
				peptide.getAccessions(),
				-1,
				parameters.getAAConstants().getChargedMass(peptide.getPeptideModSeq(), peptide.getPrecursorCharge()),
				peptide.getPrecursorCharge(),
				peptide.getPeptideModSeq(),
				0,
				warpedRTInSec,
				peptide.getScore(),
				masses,
				intensities,
				correlations,
				quantifiedIons,
				entry.getIonMobility(),
				parameters.getAAConstants(),
				true // force preserving peaks with non-positive intensity (like any we had to add)
		);
	}

	/**
	 * Quantify peptides for one or more files, based on previously-computed alignment and transition refinement.
	 *
	 * @param jobs one or more jobs that should be quantified
	 * @param elibFile the location where the results should be written in ELIB format (quantitative)
	 * @param alignmentElib the location from which previously-computed "alignment-only" results should be read
	 * @param parameters the parameters to use for quant (should match those used for the initial alignment exactly!)
	 */
	public static void convertElibQuantOnly(ProgressIndicator progress, List<? extends SearchJobData> jobs, File elibFile, File alignmentElib, SearchParameters parameters) {
		final Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides;
		final PeakLocationInferrerInterface inferrer;
		final ArrayList<? extends ProteinGroupInterface> proteins;
		final TObjectFloatMap<String> ticMap = new TObjectFloatHashMap<>();

		try {
			try {
				FileLogRecorder logRecorder = new FileLogRecorder(new File(elibFile.getAbsolutePath() + EncyclopediaJobData.LOG_FILE_SUFFIX));
				Logger.addRecorder(logRecorder);
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				Logger.errorLine("Error recording logs to results folder!");
				Logger.errorException(e);
			}

			try {
				final LibraryFile alignmentFile = new LibraryFile();
				try {
					alignmentFile.openFile(alignmentElib);

					passingPeptides = readPassingPeptides(alignmentFile, parameters);
					proteins = alignmentFile.getProteinGroups();

					inferrer = readInferrer(alignmentFile, jobs, parameters);

					// Read TICs from alignment ELIB; we can't just read the provided set of jobs, we want
					// all the jobs used for creation of the alignment results, so we match the normalization
					// as if they were all quantified together.
					try (
							Connection c = alignmentFile.getConnection();
							PreparedStatement ps = c.prepareStatement("SELECT Key, Value FROM Metadata WHERE Key LIKE (? || '%');")
					) {
						ps.setString(1, LibraryFile.SOURCEFILE_TIC_PREFIX);

						try (ResultSet rs = ps.executeQuery()) {
							while (rs.next()) {
								ticMap.put(
										rs.getString(1).substring(LibraryFile.SOURCEFILE_TIC_PREFIX.length()),
										Float.parseFloat(rs.getString(2))
								);
							}
						}
					}
					Logger.logLine("Found TIC values to normalize across " + ticMap.size() + " files");

				} finally {
					alignmentFile.close();
				}
			} catch (IOException | SQLException e) {
				throw new EncyclopediaException("Unable to read alignment results", e);
			}

			convertElibQuantOnly(progress, jobs, elibFile, passingPeptides, inferrer, proteins, ticMap, parameters, parameters.isIntegratePrecursors());
		} finally {
			Logger.close();
		}
	}

	/**
	 * Read the set of passing peptides from an "alignment-only" ELIB (ALIB).
	 *
	 * @param alignmentFile an open ALIB library
	 */
	static Pair<ArrayList<PercolatorPeptide>, Float> readPassingPeptides(LibraryFile alignmentFile, SearchParameters parameters) throws IOException, SQLException {
		final ArrayList<PercolatorPeptide> passingPeptides = Lists.newArrayList();
		float pi0 = Float.parseFloat(alignmentFile.getMetadata().get("pi0"));

		final String query = "SELECT" +
				" e.sourcefile, e.rtinseconds, max(p.isdecoy), e.peptidemodseq, e.precursorcharge," +
				" group_concat(p.proteinaccession, ';')," +
				" s.qvalue, s.posteriorerrorprobability" +
				" FROM entries e" +
				" JOIN peptidescores s USING (peptidemodseq, precursorcharge)" + // assume single rows per (modseq, z)
				" JOIN peptidetoprotein p USING (peptideseq)" +
				" GROUP BY peptidemodseq, precursorcharge;";

		try (Connection c = alignmentFile.getConnection()) {
			try (PreparedStatement ps = c.prepareStatement(query)) {
				try (ResultSet rs = ps.executeQuery()) {
					while (rs.next()) {
						passingPeptides.add(new PercolatorPeptide(
								PercolatorPeptide.getPSMID(
										rs.getString(1), // sourcefile
										rs.getFloat(2), // rt
										Optional.empty(),
										rs.getBoolean(3), // decoy
										rs.getString(4), // peptidemodseq
										rs.getByte(5) // charge
								),
								rs.getString(6), // proteinids
								rs.getFloat(7),  // qvalue
								rs.getFloat(8),  // PEP
								parameters.getAAConstants()
						));
					}
				}
			}
		}

		return new Pair<>(passingPeptides, pi0);
	}

	/**
	 * Read the RT alignment and transition refinement results from an "alignment-only" ELIB (ALIB).
	 *
	 * @param alignmentFile an open ALIB library
	 * @param jobs the set of files to which the RT alignment results shoudl be limited. May be null or empty to indicate "all".
	 *
	 * @return an appropriate {@code PeakLocationInferrerInterface} that can be used for quantification of some or all
	 *         jobs recorded in the ALIB.
	 *
	 * @see #getRawAlignmentPoints(SearchJobData, PeakLocationInferrerInterface) used in creating the data read by this method
	 */
	static PeakLocationInferrerInterface readInferrer(LibraryFile alignmentFile, List<? extends SearchJobData> jobs, SearchParameters parameters) throws IOException, SQLException {
		final HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentMap = Maps.newHashMap();
		final HashMap<SearchJobData, List<AlignmentDataPoint>> alignmentDataMap = Maps.newHashMap();
		final HashMap<String, Float> alignedRTInMinBySequenceMap = Maps.newHashMap();
		final HashMap<String, double[]> bestIons = Maps.newHashMap();

		try (Connection c = alignmentFile.getConnection()) {
			try (PreparedStatement ps = c.prepareStatement(
					"SELECT Library, Actual, Predicted, Delta, Probability, Decoy, PeptideModSeq" +
							" FROM retentiontimes WHERE sourcefile = ?;"
			)) {
				// Read alignment for each job
				for (SearchJobData job : jobs) {
					final List<AlignmentDataPoint> alignmentData = Lists.newArrayList();

					ps.setString(1, job.getOriginalDiaFileName());
					try (ResultSet rs = ps.executeQuery()) {
						while (rs.next()) {
							// Must be in _minutes_ to match AlternatePeakLocationInferrer; the values have
							// been converted to seconds when written to the `retentiontimes` table.
							alignmentData.add(AlignmentDataPoint.of(
									rs.getFloat(1) / 60f, // lib
									rs.getFloat(2) / 60f, // actual
									rs.getFloat(3) / 60f, // predicted
									rs.getFloat(4) / 60f, // delta
									rs.getFloat(5), // prob
									rs.getBoolean(6), // decoy
									rs.getString(7) // modseq
							));
						}
					}

					Logger.logLine(job.getOriginalDiaFileName() + " alignment points: read " + alignmentData.size() );

					if (alignmentData.stream().noneMatch(p -> Objects.nonNull(p.getPeptideModSeq()))) {
						// "Seed" experiment won't have alignment points and shouldn't be populated in the map.
						Logger.logLine("Assuming job " + job.getOriginalDiaFileName() + " is the seed; using 1-1 RT mapping.");
						continue;
					}

					alignmentDataMap.put(job, alignmentData);

					// Generate the necessary linear interpolation for prediction -- MUST use ALL points
					final ArrayList<XYPoint> alignmentPoints = alignmentData.stream()
							// RTs must be in _minutes_ to match AlternatePeakLocationInferrer, but these values are
							// already in minutes when we read them from `retentiontimes`.
							.map(p -> new XYPoint(p.getLibrary(), p.getPredictedActual()))
							// We sort by increasing X, but then by increasing Y, to avoid causing issues
							// with the monotonicity check below when the same X value is repeated (which
							// we see happen sometimes).
							.sorted(Comparator.comparingDouble(XYPoint::getX)
									          .thenComparing(XYPoint::getY)
							)
							.collect(Collectors.toCollection(ArrayList::new));

					// Check for monotonic function (don't require strict monotonicity though, as we see this sometimes).
					// This non-strictness means the function may not be correctly invertible! This isn't our problem
					// though, we're just handling the behavior that's in place elsewhere. We do check however to be sure
					// that the provided data isn't entirely nonsensible.
					for (int i = 1; i < alignmentPoints.size(); i++) {
						final double y = alignmentPoints.get(i).getY();
						final double prev = alignmentPoints.get(i - 1).getY();

						if (y < prev){
							throw new IllegalStateException(String.format(
									"Alignment warp is not monotonic! (%.02f, %.02f) -> (%.02f, %.02f)",
									alignmentPoints.get(i - 1).getX(),
									prev,
									alignmentPoints.get(i).getX(),
									y
							));
						}
					}

					final LinearInterpolatedFunction rtWarper = new LinearInterpolatedFunction(alignmentPoints);

					/**
					 * Points are sorted by delta (x), but function is not monotonic. NOT INVERTIBLE!
					 */
					final LinearInterpolatedFunction probFn = new LinearInterpolatedFunction(alignmentData.stream()
							// Filter out raw knots without delta/probability fields
							.filter(p -> Float.isFinite(p.getDelta()) && Floats.isFinite(p.getProbability()))
							.map(p -> new XYPoint(p.getDelta(), p.getProbability()))
							.sorted(Comparator.comparingDouble(XYPoint::getX))
							.collect(Collectors.toCollection(ArrayList::new))
					);

					Optional<RTProbabilityModel> probModel = Optional.of(new RTProbabilityModel() {
						@Override
						public float getProbability(float retentionTime, float delta) {
							return probFn.getYValue(delta);
						}
					});

					final RetentionTimeAlignmentInterface alignment = new AbstractRetentionTimeFilter(rtWarper, probModel, null, null) {
						static final double MATCH_TOLERANCE = 1e-3;

						@Override
						public List<AlignmentDataPoint> plot(List<XYPoint> rts, Optional<File> saveFileSeed) {
							return alignmentData.stream()
									.filter(p -> Objects.nonNull(p.getPeptideModSeq()))
									.collect(Collectors.toList());
						}

						@Override
						public float getYValue(float xrt) {
							return rtWarper.getYValue(xrt);
						}

						@Override
						public float getXValue(float yrt) {
							return rtWarper.getXValue(yrt);
						}

						@Override
						public float getDelta(float actualRT, float modelRT) {
							// quick-n-dirty -- require exact match
							for (AlignmentDataPoint p : alignmentData) {
								if (
										Floats.isFinite(p.getActual())
												&& Math.abs(actualRT - p.getActual()) < MATCH_TOLERANCE
												&& Math.abs(modelRT - p.getPredictedActual()) < MATCH_TOLERANCE
								) {
									return p.getDelta();
								}
							}

							throw new IllegalStateException("No such alignment point with actual/model RT " + actualRT + " / " + modelRT);
						}
					};

					alignmentMap.put(job, alignment);
				}
			}

			// Read RT, ions from entries table
			try (PreparedStatement ps = c.prepareStatement(
					"SELECT" +
					" e.PeptideModSeq," +
					" e.RTInSeconds," +
					" e.MassArray," +
					" e.MassEncodedLength," +
					" e.QuantifiedIonsArray" +
					" FROM entries e;"
			)) {
				try (ResultSet rs = ps.executeQuery()) {
					while (rs.next()) {
						final String modSeq = rs.getString(1);

						final double[] masses = ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(3), rs.getInt(4)));

						final boolean[] quantifiedIons;
						final byte[] compressedQuantIons = rs.getBytes(5);
						if (null == compressedQuantIons || 0 == compressedQuantIons.length) {
							quantifiedIons = new boolean[masses.length];
							Arrays.fill(quantifiedIons, true); // old results shouldn't be used as alignments, but OK
						} else {
							quantifiedIons = ByteConverter.toBooleanArray(CompressionUtils.decompress(compressedQuantIons));
						}

						final TDoubleList quantIonMasses = new TDoubleArrayList();
						for (int i = 0; i < masses.length; i++) {
							if (i < quantifiedIons.length) { // clamp in bounds; should be unnecessary
								if (quantifiedIons[i]) {
									quantIonMasses.add(masses[i]);
								}
							}
						}

						bestIons.put(modSeq, quantIonMasses.toArray());

						final float rtInSec = rs.getFloat(2);
						if (rs.wasNull() || !Float.isFinite(rtInSec) || rtInSec < 0) { // null/nan, infinite, or negative values should not be recorded in the alignment
							continue;
						}

						alignedRTInMinBySequenceMap.put(modSeq, rtInSec / 60f); // must be converted to minutes for use by the inferrer
					}
				} catch (DataFormatException e) {
					throw new EncyclopediaException("Invalid mass encoding!", e);
				}
			}
		}

		return new SimplePeakLocationInferrer(
				alignmentMap,
				alignmentDataMap,
				alignedRTInMinBySequenceMap,
				bestIons,
				parameters
		);
	}

	/**
	 * Quantify peptides for one or more files, based on previously-computed alignment and transition refinement.
	 *
	 * TODO: eventually this method should be combined with {@link #convertElib(ProgressIndicator, List, File, Optional, Optional, Optional, SearchParameters)},
	 *       likely by modifying that method to call this one
	 *
	 * @param jobs one or more jobs that should be quantified
	 * @param elibFile the location where the results should be written in ELIB format (quantitative)
	 * @param passingPeptides the previously-computed set of passing peptides
	 * @param inferrer the previously-computed RT alignment and transition refinement
	 * @param proteins the previously-computed set of scored and grouped proteins
	 * @param parameters the parameters to use for quant (should match those used for the initial alignment exactly!)
	 */
	static void convertElibQuantOnly(ProgressIndicator progress, List<? extends SearchJobData> jobs, File elibFile, Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides, PeakLocationInferrerInterface inferrer, ArrayList<? extends ProteinGroupInterface> proteins, TObjectFloatMap<? super String> ticMap, SearchParameters parameters, boolean integratePrecursors) {
		if (null == passingPeptides || null == passingPeptides.x || passingPeptides.x.size() < 1) {
			throw new IllegalArgumentException("Can't extract quantities for zero peptides!");
		}

		Objects.requireNonNull(proteins, "Unable to proceed without previously-computed protein groups!");

		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();

			float increment=1.0f/jobs.size();
			for (int i=0; i<jobs.size(); i++) {
				SearchJobData job = jobs.get(i);
				if (!job.hasBeenRun()) {
					Logger.errorLine("Unable to process " + job.getOriginalDiaFileName() + " because its results are missing. Continuing.");
					continue;
				}
				ProgressIndicator subProgress = new SubProgressIndicator(progress, increment);

				ArrayList<PercolatorPeptide> globalPassingPeptides = passingPeptides.x;
				Pair<ArrayList<PercolatorPeptide>, Float> localPassingPeptides = PercolatorReader.getPassingPeptidesFromTSV(job.getPercolatorFiles().getPeptideOutputFile(), jobs.get(i).getParameters(), false);

				Logger.logLine(job.getOriginalDiaFileName() + ": Number of global peptides: " + globalPassingPeptides.size() + " vs local peptides: " + localPassingPeptides.x.size());

				convertFileElib(subProgress, job, globalPassingPeptides, localPassingPeptides.x, Optional.of(inferrer), elib, jobs.size() > 1, integratePrecursors);
			}

			writeElibMetadata(elib, jobs, parameters, true);

			// Write _all_ passing peptide scores to results ELIB, as they're used downstream (quant reports)
			// without caring which sample each peptide was scored in. Don't bother with decoys (should be fine).
			elib.addTargetDecoyPeptides(passingPeptides.x, Lists.newArrayList());

			// We must also take the protein connections from the entries so they're available to quant reporting.
			elib.addProteinsFromPercolator(passingPeptides.x);

			elib.createIndices();
			elib.saveAsFile(elibFile);

			try {
				IntensityNormalizer tic = parameters.isNormalizeByTIC()?IntensityNormalizer.tic(ticMap):IntensityNormalizer::identity;
				LibraryReportExtractor.extractMatrix(
						elib,
						jobs.stream().map(j -> j.getOriginalDiaFileName()).collect(Collectors.toList()),
						proteins,
						tic,
						Optional.empty(),
						""
				);
			} catch (DataFormatException e) {
				Logger.errorException(e);
			}

			elib.close();
		} catch (IOException ioe) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error creating BLIB file");
			Logger.errorException(sqle);
		}
	}

	private static void logScribeTwoPerf(String message) {
		if (!ScribeTwo.SCRIBE2_PERF_LOGGING) return;
		Logger.logLine("SCRIBE2_PERF phase="+message+" heapUsedMiB="+getUsedHeapMiB());
	}

	private static long getUsedHeapMiB() {
		Runtime runtime=Runtime.getRuntime();
		return Math.round((runtime.totalMemory()-runtime.freeMemory())/1024.0/1024.0);
	}

}
