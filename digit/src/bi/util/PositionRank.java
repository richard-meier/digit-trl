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

public class PositionRank {
	private int position;
	private PositionRank previous, next;
	private String chromosomeName;
	private ReadTag tag;
	private PositionRank mate;
	
	public PositionRank(int position, ReadTag tag, String chromosomeNumber){
		this.position=position;
		this.setPrevious(null);
		this.setNext(null);
		this.chromosomeName=chromosomeNumber;
		this.tag=tag;
		this.setMate(null);
	}

	public PositionRank getPrevious() {
		return previous;
	}

	public void setPrevious(PositionRank previous) {
		this.previous = previous;
	}

	public PositionRank getNext() {
		return next;
	}

	public void setNext(PositionRank next) {
		this.next = next;
	}
	
	public int getPosition(){
		return position;
	}
	
	public ReadTag getTag(){
		return tag;
	}
	
	private static int testA=0, testB=0;
	public ReadTag getTagAndAddPositionalInformationToIt(){
		if(tag.getPositionTag()==null){
			tag.setPositionTag(new PositionTag(chromosomeName,position));
		}
		else{
		}
		if(this instanceof PositionRankSD){
			tag.setStrand(((PositionRankSD)(this)).getStrand());
		}
		return tag;
	}

	public String getChromosomeName() {
		return chromosomeName;
	}

	public PositionRank getMate() {
		return mate;
	}

	public void setMate(PositionRank mate) {
		this.mate = mate;
	}
}
