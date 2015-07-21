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

public class DnaMod {
	public static String reverseComplement(String dna){
		StringBuilder out = new StringBuilder();
		for(int i=dna.length()-1; i>=0; i--){
			out.append(getComplement(dna.charAt(i)));
		} 
		return out.toString();
	}
	
	private static char getComplement(char c){
		if(c == 'A') return 'T';
		if(c == 'T') return 'A';
		if(c == 'C') return 'G';
		if(c == 'G') return 'C';
		if(c == 'N') return 'N';
		else {
			System.err.println("Error:DnaMod:getComplement: Unexpected character "+c);
			System.exit(-1);
			return ' '; 
		}
	}
}
