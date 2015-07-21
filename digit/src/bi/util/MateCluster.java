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

import java.util.Vector;

public class MateCluster {
	private int[] region;
	private MateCluster mate;
	private Vector<ReadTag> occupants;
	private int concordantCounter, discordantCounter;
	private String chromosomeName;
	private int[] arrayForCoverageCalculation;
	
	public MateCluster(int start, int stop, String chrName){
		region = new int[2];
		region[0]=start;
		region[1]=stop;
		this.chromosomeName=chrName;
		mate=null;
	}
	
	public int getStart(){return region[0];}
	public int getStop(){return region[1];}
	public MateCluster getMate(){return mate;}
	
	public void setMate(MateCluster mate){this.mate=mate;}

	public Vector<ReadTag> getOccupants() {
		return occupants;
	}
	
	public void initilizeOccupants() {
		occupants = new Vector<ReadTag>();
	}

	@SuppressWarnings("unchecked")
	public void addOccupants(Vector<ReadTag> occupants) {
		if(occupants == null) this.occupants=null;
		else this.occupants = (Vector<ReadTag>)(occupants.clone());
	}
	
	public void addOccupant(ReadTag read) {
		occupants.add(read);
	}
	
	public String getOccupantsAsString(){
		String line = "";
		for(int i=0;i<occupants.size();i++){
			line+=occupants.get(i).getName()+occupants.get(i).getFileNumber()+'\t';
		}
		return line;
	}
	
	public void appendDescription(StringBuffer write, boolean appendReadNames){
		write.append(hashCode()); write.append("\t");
		write.append(chromosomeName); write.append("\t");
		write.append(region[0]); write.append("\t");
		write.append(region[1]); write.append("\t");
		write.append(mate.hashCode()); write.append("\t");
		if(appendReadNames) appendOccupantNames(write);
		else write.append("-");
	}
	
	public void appendCircosEntryTo(StringBuffer write){
		write.append(chromosomeName); write.append(" ");
		write.append(region[0]); write.append(" ");
		write.append(region[1]); write.append(" ");
	}
	
	public void appendOccupantNames(StringBuffer write){
		for(int i=0; i<occupants.size()-1; i++) {
			ReadTag t=occupants.get(i);
			t.appendFullName(write);
			write.append("\t");
		}
		occupants.get(occupants.size()-1).appendFullName(write);
	}

	public int getConcordantCounter() {
		return concordantCounter;
	}

	public void setConcordantCounter(int concordantCounter) {
		this.concordantCounter = concordantCounter;
	}

	public void increaseConcordantCounter() {
		concordantCounter++;
	}
	
	public int getDiscordantCounter() {
		return discordantCounter;
	}

	public void setDiscordantCounter(int discordantCounter) {
		this.discordantCounter = discordantCounter;
	}

	public String getChromosomeName() {
		return chromosomeName;
	}

	public void setChromosomeName(String chromosomeName) {
		this.chromosomeName = chromosomeName;
	}
	
	public boolean doesOverlapWith(MateCluster otherCluster){
		return (otherCluster.getStop()>=this.getStart() && otherCluster.getStart()<=this.getStart()) ||
			(this.getStop()>=otherCluster.getStop() && this.getStart()<=otherCluster.getStart())
		;
	}
	
	public boolean containsAsAWhole(MateCluster otherCluster){
		if(otherCluster.getChromosomeName().equals(this.getChromosomeName())){
			if(otherCluster.getStart()>=this.getStart()){
				if(otherCluster.getStop()<=this.getStop()) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	}
	
	public boolean containsAsAWhole(int pos, String chr){
		if(chr.equals(this.getChromosomeName())){
			if(pos >= this.getStart()){
				if(pos <= this.getStop()) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	}
	
	public boolean isInClusterRange(int pos, int outerRange, String chr){
		if(chr.equals(this.getChromosomeName())){
			if(pos >= this.getStart()-outerRange){
				if(pos <= this.getStop()+outerRange) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	}
	
	public void merge(MateCluster toBeIncluded){
		for(ReadTag t:toBeIncluded.getOccupants()){
			if(!occupants.contains(t)) addOccupant(t);
		}
		if(toBeIncluded.getStart()<region[0]) region[0]=toBeIncluded.getStart();
		if(toBeIncluded.getStop()>region[1]) region[1]=toBeIncluded.getStop();
	}

	public void inizializeArrayForCoverageCalculation(int size) {
		arrayForCoverageCalculation = new int [size];
	}
	
	public void addReadToArrayForCoverageCalculation(int start, int length) {
		if(start<0){
			length+=start;
			start=0;
		}
		for(int i=0; i<length;i++){
			if((i+start)>=arrayForCoverageCalculation.length)break;
			arrayForCoverageCalculation[i+start]++;
		}
	}
	
	public int getClusterSize() {
		return getStop()-getStart();
	}

	public double getCoverage() {
		int sum=0;
		for(int i=0;i<arrayForCoverageCalculation.length;i++)
			sum+=arrayForCoverageCalculation[i];
		return (1.0*sum)/arrayForCoverageCalculation.length;
	}
}
