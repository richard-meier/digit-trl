/**
 * Copyright 2014 Richard Meier, Stefan Graw, Jeremy Chien, 
 * Peter Beyerlein
 *
 *  This file is part of the software pipeline digit.
 *
 *  digit is free software: you can redistribute it and/or modify it 
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  digit is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *  
 *  For a copy of the GNU General Public License see 
 *  <http://www.gnu.org/licenses/>.
 *
 */

package bi.util;


public class Cluster {

	private MateCluster mate1, mate2;
	private String chr1, chr2;
	 
	
	public Cluster(String input){
		String [] inputArray = input.split("\\s");
		setChr1(hsToChr(inputArray[0]));
		int startR1 = Integer.parseInt(inputArray[1]);
		int stopR1 = Integer.parseInt(inputArray[2]);	
		setChr2(hsToChr(inputArray[3]));
		int startR2 = Integer.parseInt(inputArray[4]);
		int stopR2 = Integer.parseInt(inputArray[5]);
		
		setMate1(new MateCluster(startR1, stopR1, chr1));
		setMate2(new MateCluster(startR2, stopR2, chr2));
		mate1.setMate(mate2);
		mate2.setMate(mate1);
	}

	public Cluster(MateCluster mate1, MateCluster mate2){
		this.mate1=mate1;
		this.mate2=mate2;
		chr1=mate1.getChromosomeName();
		chr2=mate2.getChromosomeName();
	}

	protected String hsToChr(String hs) {
		String chr = "chr";
		for(int i=2; i<hs.length();i++)
			chr+=hs.charAt(i);
		return chr;
	}


	public String getChr1() {
		return chr1;
	}


	public void setChr1(String chr1) {
		this.chr1 = chr1;
	}


	public String getChr2() {
		return chr2;
	}


	public void setChr2(String chr2) {
		this.chr2 = chr2;
	}


	public MateCluster getMate1() {
		return mate1;
	}


	public void setMate1(MateCluster mate1) {
		this.mate1 = mate1;
	}

	public MateCluster getMate2() {
		return mate2;
	}

	public void setMate2(MateCluster mate2) {
		this.mate2 = mate2;
	}
	
	public boolean contains(MateCluster mc){
		return mate1==mc || mate2==mc;
	}

	public boolean containsSubdivision(Cluster subdivision) {
		MateCluster devisionMate1 = subdivision.getMate1();
		MateCluster devisionMate2 = subdivision.getMate2();
		if( mate1.containsAsAWhole(devisionMate1) && mate2.containsAsAWhole(devisionMate2) ){
			return true;
		}
		else if( mate1.containsAsAWhole(devisionMate2) && mate2.containsAsAWhole(devisionMate1) ){
			return true;
		}
		else return false;
	}
	
	public boolean containsPosition(int pos, String chr){
		return mate1.containsAsAWhole(pos, chr) || mate2.containsAsAWhole(pos, chr);
	}
	
	public MateCluster isInClusterRange(int pos, int outerRange, String chr){
		if(mate1.isInClusterRange(pos, outerRange, chr)){
			return mate1;
		}
		else if(mate2.isInClusterRange(pos, outerRange, chr)){
			return mate2;
		}
		else return null;
	}
	
	public String toString(){
		StringBuffer out = new StringBuffer();
		mate1.appendCircosEntryTo(out);
		mate2.appendCircosEntryTo(out);
		out.append("\n");
		return out.toString();
	}
	
	public int getLengthSum(){
		return mate1.getStop()-mate1.getStart()+mate2.getStop()-mate2.getStart();
	}
	
	public MateCluster getMateOnChromosome(String chr){
		if(mate1.getChromosomeName().equals(chr)){
			return mate1;
		}
		else{
			return mate2;
		}
	}
}
