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

public class SAMEntry implements Comparable<SAMEntry>{

	private String qname, rname, quality, sequence;
	private int flag, position;
	private double mappingQuality;
	private String tagXT,tagXS,tagZM,cigar;
	
	public SAMEntry(String line){
		String[] entries = line.split("\\s");
		setQname(entries[0]);
		setRname(entries[2]);
		setFlag(Integer.parseInt(entries[1]));
		setPosition(Integer.parseInt(entries[3]));
		setMappingQuality(Double.parseDouble(entries[4]));
		setCigar(entries[5]);
		setSequence(entries[9]);
		setQuality(entries[10]);
		tagXT=null;
		for(String e:entries){
			if(e.contains("XT:A:")){
				tagXT=e;
				continue;
			}
			if(e.contains("XS:i:")){
				tagXS=e;
				continue;
			}
			if(e.contains("ZM:Z:")){
				tagZM=e;
				continue;
			}
		}
	}

	public String getZM() {
		return tagZM;
	}
	
	public String getXT() {
		return tagXT;
	}
	
	public String getXS() {
		return tagXS;
	}
	
	public String getQname() {
		return qname;
	}

	public void setQname(String qname) {
		this.qname = qname;
	}

	public String getRname() {
		return rname;
	}

	public void setRname(String rname) {
		this.rname = rname;
	}

	public int getFlag() {
		return flag;
	}

	public void setFlag(int flag) {
		this.flag = flag;
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public String getQuality() {
		return quality;
	}

	public void setQuality(String quality) {
		this.quality = quality;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public double getMappingQuality() {
		return mappingQuality;
	}

	public void setMappingQuality(double mappingQuality) {
		this.mappingQuality = mappingQuality;
	}

	public String getCigar() {
		return cigar;
	}

	public void setCigar(String cigar) {
		this.cigar = cigar;
	}

	@Override
	public int compareTo(SAMEntry o) {
		if(o.getPosition()<this.getPosition()){
			return -1;
		}
		else if(o.getPosition()>this.getPosition()){
			return 1;
		}
		else{
			return 0;
		}
	}
}
