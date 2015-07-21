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

public class AlignmentResult {
	public String alignmentString1,alignmentString2,startString1,startString2;
	public double score;
	public int sequenceStart1, sequenceStart2, sequenceEnd1, sequenceEnd2;
	
	public String toString(){
		StringBuffer out = new StringBuffer();
		out.append(alignmentString1);
		out.append("\n");
		out.append(alignmentString2);
		out.append("\n");
		out.append(sequenceStart1);
		out.append("\t");
		out.append(sequenceEnd1);
		out.append("\n");
		out.append(sequenceStart2);
		out.append("\t");
		out.append(sequenceEnd2);
		out.append("\n");
		out.append(score);
		out.append("\n");
		return out.toString();
	}
	
	public double lNormScore1(){
		return score * (sequenceEnd1-sequenceStart1)/(startString1.length()*1.0);
	}
	
	public double lNormScore2(){
		return score * (sequenceEnd2-sequenceStart2)/(startString2.length()*1.0);
	}
	
	public double lFrac1(){
		return (sequenceEnd1-sequenceStart1)/(startString1.length()*1.0);
	}
	
	public double lFrac2(){
		return (sequenceEnd2-sequenceStart2)/(startString2.length()*1.0);
	}
}
