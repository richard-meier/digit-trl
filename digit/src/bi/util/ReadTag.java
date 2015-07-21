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


public class ReadTag{

	private String name;
	private byte fileNumber;
	private PositionTag pt;
	private char strand;

	public ReadTag(String name, byte fileNumber){
		this.name=name;
		this.fileNumber=fileNumber;
		this.strand='|';
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public byte getFileNumber() {
		return fileNumber;
	}
	public void setFileNumber(byte fileNumber) {
		this.fileNumber = fileNumber;
	}
	
	public void setPositionTag(PositionTag pt){
		this.pt=pt;
	}
	
	public PositionTag getPositionTag(){
		return pt;
	}
	
	public char getStrand(){
		return strand;
	}
	
	public boolean equals(ReadTag other){
		if(other.getName().equals(name) && other.getFileNumber()==fileNumber) return true;
		else return false;
	}
	
	public void appendFullName(StringBuffer write){
		write.append(name);
		write.append(fileNumber);
		write.append(":"+strand+":");
		if(pt!=null){
			write.append(":"+pt.getPosition());
		}
	}

	public void setStrand(char strand) {
		this.strand=strand;
	}
}
