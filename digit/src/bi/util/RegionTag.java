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

public class RegionTag {
	private String chromosome;
	private int start, stop;
	
	public RegionTag(String chromosome, int start, int stop){
		this.setChromosome(chromosome);
		this.setStart(start);
		this.setStop(stop);
	}
	
	public RegionTag(String[] information){
		this(information[0],Integer.parseInt(information[1]),Integer.parseInt(information[2]));
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getStop() {
		return stop;
	}

	public void setStop(int stop) {
		this.stop = stop;
	}
	
	public int distanceFromPoint(int i){
		if(i<start){
			return i-start;
		}
		else if(i>stop){
			return i-stop;
		}
		else return 0;
	}
}
