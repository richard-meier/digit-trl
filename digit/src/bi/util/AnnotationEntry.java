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

public class AnnotationEntry {

	private String geneID, geneName, chromosomeID, feature;
	private int start,stop;
	private char strand;
	private String type;
	
	public AnnotationEntry(String line){
		String[] columns = line.split("\\t");
		chromosomeID=columns[0];
		feature=columns[2];
		start=Integer.parseInt(columns[3]);
		stop=Integer.parseInt(columns[4]);
		strand=columns[6].charAt(0);
		setType(columns[1]);
		String attributeString = columns[8];
		attributeString = attributeString.replaceAll(";{1}\\s{1}", ";");
		attributeString = attributeString.replaceAll("\\s", "=");
		attributeString = attributeString.replaceAll("\"","");
		String[] attributes = attributeString.split(";");
		for(String s:attributes){
			String[] slots = s.split("=");
			if(slots[0].equals("gene_id")){
				geneID=slots[1];
			}
			if(slots[0].equals("gene_name")){
				setGeneName(slots[1]);
			}
		}
	}
	
	public void setGeneID(String newID) {
		geneID=newID;
	}
	
	public String getGeneID() {
		return geneID;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public String getChromosomeID() {
		return chromosomeID;
	}

	public String getFeature() {
		return feature;
	}
	
	public char getStrand(){
		return strand;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}
}

