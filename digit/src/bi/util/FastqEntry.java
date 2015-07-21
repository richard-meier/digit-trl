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

public class FastqEntry {
	private String head, sequence, description, quality;

	public FastqEntry(String head, String sequence, String description, String quality){
		this.setHead(head);
		this.setSequence(sequence);
		this.setDescription(description);
		this.setQuality(quality);
	}
	
	public FastqEntry(String fastq){
		fastq.split("");
	}

	public String getHead() {
		return head;
	}

	public String getSequence() {
		return sequence;
	}

	public String getDescription() {
		return description;
	}

	public String getQuality() {
		return quality;
	}

	public void setHead(String head) {
		this.head = head;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public void setQuality(String quality) {
		this.quality = quality;
	}
	
	public String toString(){
		return head+"\t"+sequence+"\t"+description+"\t"+quality+"\t";
	}
	
	public String toFileFormatString(){
		return head+"\n"+sequence+"\n"+description+"\n"+quality+"\n";
	}
}
