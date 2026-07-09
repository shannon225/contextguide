package org.searlelab.contextguide.mprophet;

public class IsolationWindow {
	private double targetMz;
	private double mzStart;
	private double mzStop; 
	private double windowMz;
	private float rtCenter; // Including this because we use it to calculate m/z min and max in the TargetedBoostrapper
	private float rtMin;
	private float rtMax;
	private boolean isDecoy;
	private String compound;
	private byte charge;
	private String peptide;

	
	// Constructor
	public IsolationWindow(double precursorMz, float rtMin, float rtMax, boolean isDecoy) {
		this.targetMz = precursorMz;
		this.rtMin = rtMin;
		this.rtMax = rtMax;
		this.isDecoy = isDecoy;
		}
	
	// Isolation window as mass list with peptide and charge state
	public IsolationWindow(String peptide, double precursorMz, byte charge, float rtMin, float rtMax, boolean isDecoy) {
		this.peptide = peptide;
		this.targetMz = precursorMz;
		this.charge = charge;
		this.rtMin = rtMin;
		this.rtMax = rtMax;
		this.isDecoy = isDecoy;
	}
	
	// Getters 
	public double getTargetMz() {
		return targetMz;
	}
	
	public double getWindowMz() {
		return windowMz;
	}
	
	public double getMzStart() {
		return mzStart;
	}
	
	public double getMzStop() {
		return mzStop;
	}

	public float getRtCenter() {
		return rtCenter;
	}
	
	public float getRtMin() {
		return rtMin;
	}
	
	public float getRtMax() {
		return rtMax;
	}
	
	public boolean isDecoy() {
		return isDecoy;
	}

	public String getCompound() {
		return compound;
	}

	public byte getCharge() {
		return charge;
	}
	
	public String getPeptide() {
		return peptide;
	}

}

