package org.searlelab.context.encyclopedia;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetDataset;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.ScoredMProphetData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class MProphetReiter implements Runnable {
	private final float peptideFDRThreshold;
	private final MProphetExecutionData settings;
	private final AminoAcidConstants aaConstants;
	private final int seed;
	
	private Throwable error;
	private MProphetResult result;
	
	public MProphetReiter(MProphetExecutionData settings, float peptideFDRThreshold, int seed, AminoAcidConstants aaConstants) {
		this.settings = settings;
		this.peptideFDRThreshold=peptideFDRThreshold;
		this.aaConstants=aaConstants;
		this.seed=seed;
	}
	
	// Run the excuteMProphetTSV function with a fixed seed
	public static MProphetResult executeMProphetTSV(MProphetExecutionData commandData, float threshold, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		int seed=1;
		return executeMProphetTSV(commandData, threshold, seed, aaConstants, round);
	}
	
	// Run the execute mPRophet with seed as an input parameter
	public static MProphetResult executeMProphetTSV(MProphetExecutionData commandData, float threshold, int seed, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphetReiter prophet=new MProphetReiter(commandData, threshold, seed, aaConstants);
		prophet.run();
		return prophet.getResult();
	}
	
	
	public void runWithModel(LinearDiscriminantAnalysis lda) {
		File file = settings.getInputTSV();
		
		try {
			MProphetDataset data = MProphetFeatureReader.parseFeatureFile(file, settings);
			result = calculateProbabilitiesWithModel(data, lda);
		} catch (Throwable t ) {
			Logger.errorLine("Error performing MProphetReiter with the supplied LDA model");
			Logger.errorException(t);
			this.error = t;
		}
	}
	public static MProphetResult executeMProphetTSVWithModel(MProphetExecutionData commandData, float threshold, LinearDiscriminantAnalysis lda, AminoAcidConstants aaConstants) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphetReiter prophet = new MProphetReiter(commandData, threshold, 1, aaConstants);
		prophet.runWithModel(lda);
		return prophet.getResult();
	}

	@Override
	public void run() {
		File file=settings.getInputTSV();

		try {
			MProphetDataset data = MProphetFeatureReader.parseFeatureFile(file, settings);
			result = calculateProbabilities(data, seed);

		} catch (Throwable t) {
			Logger.errorLine("Error performing mProphet!");
			Logger.errorException(t);

			this.error = t;
		}
	}

	protected MProphetResult calculateProbabilities(MProphetDataset dataset, int seed) throws EncyclopediaException {
		int randomSeed=RandomGenerator.randomInt(seed);
		int iterationCount=50;
		int maxKeptModels=iterationCount/2;
		int numIterationsPerCalculation=10;
		
		TObjectDoubleHashMap<String> seedCoefficients=new TObjectDoubleHashMap<String>();
		seedCoefficients.put("HyperScore", 0.6);
		seedCoefficients.put("xCorrLib", 5.0);
		seedCoefficients.put("xCorrModel", 0.2);
		seedCoefficients.put("numberOfMatchingPeaksAboveThreshold", 0.2);
		seedCoefficients.put("isotopeDotProduct", 0.5);
		seedCoefficients.put("correlationToPrecursor", 0.4);
		seedCoefficients.put("isIntegratedSignal", 0.2);
		seedCoefficients.put("numPeaksWithGoodCorrelation", 0.1);
		
		seedCoefficients.put("main_var_Intensity", 1.0); // W0
		seedCoefficients.put("var_Shape_(weighted)", 4.0); // W4
		seedCoefficients.put("var_Co-elution_(weighted)", -0.05); // W5
		seedCoefficients.put("var_Co-elution_count", 1.0); // W1
		
		
		double[] coefficients=new double[dataset.getFeatureNames().size()];
		for (int i = 0; i < coefficients.length; i++) {
			double value=seedCoefficients.get(dataset.getFeatureNames().get(i));
			if (value!=0.0) {
				coefficients[i]=value;
			}
		}
		LinearDiscriminantAnalysis seedModel=new LinearDiscriminantAnalysis(coefficients, 0.0);
		
		ArrayList<ScoredObject<LinearDiscriminantAnalysis>> models=new ArrayList<ScoredObject<LinearDiscriminantAnalysis>>();
		for (int n = 0; n < iterationCount; n++) {
			randomSeed=RandomGenerator.randomInt(randomSeed);
			MProphetDataset[] folds=MProphetDataset.splitKFold(dataset, 2, randomSeed, settings.getParameters().getPercolatorTrainingSetSize());
			
			MProphetDataset trainingDataset=folds[0];
			MProphetDataset testingDataset=folds[1];
			
			LinearDiscriminantAnalysis lda;
			if(General.checkNaN(seedModel.getCoefficients())||General.checkIfAllX(seedModel.getCoefficients(), 0.0)) {
				lda=null;
			} else {
				lda=seedModel;
			}
			
			int best=0;
			LinearDiscriminantAnalysis previousLDA=null;
			for (int i = 0; i < numIterationsPerCalculation; i++) {
				float targetFDR=0.01f;
				if (i==0) {
					targetFDR=0.15f; // first time through Skyline starts with 15% of the target data, not a 15% FDR. This ensures that the first round of training can happen!
				} else if (i==1) {
					targetFDR=0.02f;
				}
				
				ArrayList<ScoredMProphetData> data;
				if (i==0) {
					data=trainingDataset.getPassingTargetsByPercentage(Optional.ofNullable(lda), targetFDR);
				} else {
					data=trainingDataset.getPassingTargetsByFDR(Optional.ofNullable(lda), targetFDR).x;
				}
				
				System.out.println(n+"."+i+") "+data.size()+" (targetFDR="+Math.round(targetFDR*100)+"%) starting model: "+lda); //FIXME
				if (data.size()==0||(i>2&&data.size()<=best)) {
					if (previousLDA!=null) {
						lda=previousLDA;
					}
					break;
				}
				best=data.size();
				ArrayList<float[]> decoyData = trainingDataset.getDecoyData();
				ArrayList<float[]> scoredData = MProphetDataset.getScoredData(data);
				if (decoyData.size()>scoredData.size()*10) {
					// make sure the decoy data doesn't get too out of hand in size
					decoyData=new ArrayList<float[]>(decoyData.subList(0, scoredData.size()*10));
				}
				LinearDiscriminantAnalysis model=LinearDiscriminantAnalysis.buildModel(scoredData, decoyData);
				
				if (!General.checkNaN(model.getCoefficients())&&!General.checkIfAllX(model.getCoefficients(), 0.0)) {
					// if we get NaNs, then fall back on wherever we were previously
					previousLDA=lda;
					lda=model;
				} else {
					break;
				}
			}
			
			if (lda==null) {
				Logger.logLine("Iteration "+(n+1)+": Failed to generate a meaningful model!");
			} else {
				Pair<ArrayList<ScoredMProphetData>, Float> data=testingDataset.getPassingTargetsByFDR(Optional.ofNullable(lda), 0.01f);
				Pair<ArrayList<ScoredMProphetData>, Float> seedData=testingDataset.getPassingTargetsByFDR(Optional.ofNullable(seedModel), 0.01f);
				if (seedData.x.size()>data.x.size()) {	
					models.add(new ScoredObject<LinearDiscriminantAnalysis>(seedData.x.size(), seedModel));
					Logger.logLine("Iteration "+(n+1)+": prefer seed model, "+seedData.x.size()+"/"+testingDataset.getTargetData().size()+" passing, pi0:"+seedData.y);
				} else {
					models.add(new ScoredObject<LinearDiscriminantAnalysis>(data.x.size(), lda));
					Logger.logLine("Iteration "+(n+1)+": "+data.x.size()+"/"+testingDataset.getTargetData().size()+" passing, pi0:"+data.y);
				}
			}
		}
		
		LinearDiscriminantAnalysis averageModel;
		if (models.size()==0) {
			Logger.logLine("No meaningful models generated, falling back on seed model for separation!");
			averageModel=seedModel;
		} else {	
			Collections.sort(models);
			Collections.reverse(models);
			// keep the N best models
			ArrayList<LinearDiscriminantAnalysis> bestModels=new ArrayList<LinearDiscriminantAnalysis>();
			for (ScoredObject<LinearDiscriminantAnalysis> scoredModel : models) {
				bestModels.add(scoredModel.getY());
				if (bestModels.size()>=maxKeptModels) break;
			}
			averageModel=LinearDiscriminantAnalysis.average(bestModels);
		}

		Pair<ArrayList<ScoredMProphetData>, Float> finalData=dataset.getPassingTargetsByFDR(Optional.ofNullable(averageModel), Float.MAX_VALUE);
		int passingCount=0;
		for (ScoredMProphetData data : finalData.x) {
			if (data.getFDR()<0.01) passingCount++;
		}
		Pair<ArrayList<ScoredMProphetData>, Float> finalDecoyData=dataset.getPassingTargetsByFDR(Optional.ofNullable(averageModel), Float.MAX_VALUE, true);

		// report model for logging
		Logger.logLine("Final model: "+passingCount+"/"+dataset.getTargetData().size()+" passing, pi0:"+finalData.y);
		
		for (int i = 0; i < averageModel.getCoefficients().length; i++) {
			Logger.logLine("   "+dataset.getFeatureNames().get(i)+" --> "+averageModel.getCoefficients()[i]);
		}
		Logger.logLine("   constant --> "+averageModel.getConstant());
		
		ArrayList<ScoredMProphetData> allData=new ArrayList<ScoredMProphetData>();
		allData.addAll(finalData.x);
		allData.addAll(finalDecoyData.x);
		Collections.sort(allData);
			
		HashSet<String> detectedPeptideSequences=new HashSet<String>();
		ArrayList<PercolatorPeptide> detectedPeptides=new ArrayList<>();
		try {
			PrintWriter targetWriter=new PrintWriter(settings.getPeptideOutputFile(), "UTF-8");
			PrintWriter decoyWriter=new PrintWriter(settings.getPeptideDecoyFile(), "UTF-8");
			
			targetWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			decoyWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			
			for (ScoredMProphetData scoredData : allData) {
				if (detectedPeptideSequences.contains(scoredData.getData().getSequence())) {
					// only keep the best (highest scoring) representation of the peptide
					continue;
				}
				detectedPeptideSequences.add(scoredData.getData().getSequence());
				
				float score=scoredData.getScore();
				float qValue=(float)scoredData.getFDR();
				float posteriorErrorProb=(float)scoredData.getLocalFDR();
				if (qValue<=peptideFDRThreshold&&!scoredData.getData().isDecoy()) {
					PercolatorPeptide pep=new PercolatorPeptide(scoredData.getData().getId(), scoredData.getData().getProtein(), qValue, posteriorErrorProb, aaConstants);
					detectedPeptides.add(pep);
				}
				
				if (scoredData.getData().isDecoy()) {
					decoyWriter.println(scoredData.getData().getId()+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.getData().getSequence()+".-"+"\t"+scoredData.getData().getProtein());
				} else {
					targetWriter.println(scoredData.getData().getId()+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.getData().getSequence()+".-"+"\t"+scoredData.getData().getProtein());
				}
			}
			targetWriter.println("pi_0="+finalData.y);
			decoyWriter.println("pi_0="+finalData.y);

			targetWriter.flush();
			decoyWriter.flush();
			targetWriter.close();
			decoyWriter.close();
			
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		}
		
		return new MProphetResult(detectedPeptides, averageModel, dataset.getFeatureNames(), finalData.y);
	}
	
	protected MProphetResult calculateProbabilitiesWithModel(
			MProphetDataset dataset,
			LinearDiscriminantAnalysis lda)
			throws EncyclopediaException {

		Pair<ArrayList<ScoredMProphetData>, Float> finalData =
				dataset.getPassingTargetsByFDR(Optional.ofNullable(lda), Float.MAX_VALUE);

		int passingCount = 0;
		for (ScoredMProphetData data : finalData.x) {
			if (data.getFDR() < 0.01) passingCount++;
		}

		Pair<ArrayList<ScoredMProphetData>, Float> finalDecoyData =
				dataset.getPassingTargetsByFDR(Optional.ofNullable(lda), Float.MAX_VALUE, true);

		Logger.logLine("Final supplied model: " + passingCount + "/"
				+ dataset.getTargetData().size() + " passing, pi0:" + finalData.y);

		for (int i = 0; i < lda.getCoefficients().length; i++) {
			Logger.logLine("   " + dataset.getFeatureNames().get(i) + " --> " + lda.getCoefficients()[i]);
		}
		Logger.logLine("   constant --> " + lda.getConstant());

		ArrayList<ScoredMProphetData> allData = new ArrayList<ScoredMProphetData>();
		allData.addAll(finalData.x);
		allData.addAll(finalDecoyData.x);
		Collections.sort(allData);

		HashSet<String> detectedPeptideSequences = new HashSet<String>();
		ArrayList<PercolatorPeptide> detectedPeptides = new ArrayList<>();

		try {
			PrintWriter targetWriter = new PrintWriter(settings.getPeptideOutputFile(), "UTF-8");
			PrintWriter decoyWriter = new PrintWriter(settings.getPeptideDecoyFile(), "UTF-8");

			targetWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			decoyWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");

			for (ScoredMProphetData scoredData : allData) {
				if (detectedPeptideSequences.contains(scoredData.getData().getSequence())) {
					continue;
				}

				detectedPeptideSequences.add(scoredData.getData().getSequence());

				float score = scoredData.getScore();
				float qValue = (float) scoredData.getFDR();
				float posteriorErrorProb = (float) scoredData.getLocalFDR();

				if (qValue <= peptideFDRThreshold && !scoredData.getData().isDecoy()) {
					PercolatorPeptide pep = new PercolatorPeptide(
							scoredData.getData().getId(),
							scoredData.getData().getProtein(),
							qValue,
							posteriorErrorProb,
							aaConstants
					);
					detectedPeptides.add(pep);
				}

				if (scoredData.getData().isDecoy()) {
					decoyWriter.println(scoredData.getData().getId() + "\t"
							+ score + "\t"
							+ qValue + "\t"
							+ posteriorErrorProb + "\t"
							+ "-." + scoredData.getData().getSequence() + ".-" + "\t"
							+ scoredData.getData().getProtein());
				} else {
					targetWriter.println(scoredData.getData().getId() + "\t"
							+ score + "\t"
							+ qValue + "\t"
							+ posteriorErrorProb + "\t"
							+ "-." + scoredData.getData().getSequence() + ".-" + "\t"
							+ scoredData.getData().getProtein());
				}
			}

			targetWriter.println("pi_0=" + finalData.y);
			decoyWriter.println("pi_0=" + finalData.y);

			targetWriter.flush();
			decoyWriter.flush();
			targetWriter.close();
			decoyWriter.close();

		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "
					+ settings.getPeptideOutputFile().getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "
					+ settings.getPeptideOutputFile().getAbsolutePath(), e);
		}

		return new MProphetResult(detectedPeptides, lda, dataset.getFeatureNames(), finalData.y);
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
	
	public MProphetResult getResult() {
		return result;
	}
}