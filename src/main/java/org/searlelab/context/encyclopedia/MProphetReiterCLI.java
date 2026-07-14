package org.searlelab.context.encyclopedia;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MProphetReiterCLI {
	public static void main(String[] args) throws Exception {
		if (args.length!=4) {
			Logger.errorLine("MProphetReiter requires four parameters in order:");
			Logger.logLine("  1) Input TSV");
			Logger.logLine("  2) Input FASTA");
			Logger.logLine("  3) Threshold (e.g., 0.01)");
			Logger.logLine("  4) Seed (e.g., 1)");
		}
		
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File inputTSV=new File(args[0]);
		File fastaFile=new File(args[1]);
		float threshold=Float.parseFloat(args[2]);
		int seed=Integer.parseInt(args[3]);

		Logger.logLine("Input TSV: "+inputTSV.toString());
		Logger.logLine("Input FASTA: "+fastaFile.toString());
		
		File peptideOutputFile=new File(inputTSV.toString()+".output.txt");
		File peptideDecoyFile=new File(inputTSV.toString()+".decoy.txt");
		
		MProphetExecutionData data=new MProphetExecutionData(inputTSV, fastaFile, peptideOutputFile, peptideDecoyFile, params);
		
		MProphetReiter.executeMProphetTSV(data, threshold, seed, params.getAAConstants(), 1);
	}
	
}
