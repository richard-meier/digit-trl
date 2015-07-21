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


public class ClusterInformation extends Cluster{
	private String type, clusterDiscrepencyReport;
	private RegionTag[] breakpointRange1,breakpointRange2;
	
	public ClusterInformation(MateCluster mate1, MateCluster mate2) {
		super(mate1, mate2);
	}
	public ClusterInformation(String line){
		super(line);
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getClusterDiscrepencyReport() {
		return clusterDiscrepencyReport;
	}
	public void setClusterDiscrepencyReport(String clusterDiscrepencyReport) {
		this.clusterDiscrepencyReport = clusterDiscrepencyReport;
	}
	public RegionTag[] getBreakpointRange1() {
		return breakpointRange1;
	}
	public RegionTag[] getBreakpointRange2() {
		return breakpointRange2;
	}
	
	public void setBreakpointRange(String[] information){
		information[0]=hsToChr(information[0]);
		if(information[0].equals(super.getMate1().getChromosomeName())){
			setBreakpointRange(information, true);
		}
		else{
			setBreakpointRange(information, false);
		}
	}
	
	private void setBreakpointRange(String[] information, boolean rangeOfMate1) {
		RegionTag[] breakpointRange;
		int iSize = information.length;
		int brSize = iSize/3;
		breakpointRange = new RegionTag[brSize];
		for(int i=0; i<brSize; i+=3){
			breakpointRange[i]=new RegionTag(information[i], Integer.parseInt(information[i+1]), Integer.parseInt(information[i+2]));
		}
		if(rangeOfMate1) breakpointRange1=breakpointRange;
		else breakpointRange2=breakpointRange;
	}
	
	public String getIdentificationString(){
		String out = "";
		out+=super.getMate1().getChromosomeName()+","+super.getMate1().getStart()+","+super.getMate1().getStop();
		out+="|";
		out+=super.getMate2().getChromosomeName()+","+super.getMate2().getStart()+","+super.getMate2().getStop();
		return out;
	}
}
