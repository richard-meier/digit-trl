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

import java.io.File;
import java.util.HashMap;
import java.util.Vector;

public class ReadPair{	
	public RegionTag reg1;
	public RegionTag reg2;
	public String information;
	public Vector<File> files;
	public Vector<String> iDs;
	public HashMap<String,Vector<ReadOrientation>> fileNameToReadsInRegion1;
	public HashMap<String,Vector<ReadOrientation>> fileNameToReadsInRegion2;
	
	public boolean overlaps(ReadPair other, int outerRange){
		boolean overlap1 = 
			(reg1.getStart()>= other.reg1.getStart()-outerRange && reg1.getStart()<= other.reg1.getStop()+outerRange) ||
			(reg1.getStop()>= other.reg1.getStart()-outerRange && reg1.getStop()<= other.reg1.getStop()+outerRange) ||
			(reg1.getStart()<= other.reg1.getStart()-outerRange && reg1.getStop()>= other.reg1.getStop()+outerRange) ||
			(reg1.getStart()>= other.reg1.getStart()-outerRange && reg1.getStop()<= other.reg1.getStop()+outerRange);
		boolean overlap2 = 
			(reg2.getStart()>= other.reg2.getStart()-outerRange && reg2.getStart()<= other.reg2.getStop()+outerRange) ||
			(reg2.getStop()>= other.reg2.getStart()-outerRange && reg2.getStop()<= other.reg2.getStop()+outerRange) ||
			(reg2.getStart()<= other.reg2.getStart()-outerRange && reg2.getStop()>= other.reg2.getStop()+outerRange) ||
			(reg2.getStart()>= other.reg2.getStart()-outerRange && reg2.getStop()<= other.reg2.getStop()+outerRange);
		return overlap1 && overlap2;
	}
	
	public boolean halfOverlaps(String chr, int otherStart, int otherStop, int outerRange){
		boolean overlap1 = chr.equals(reg1.getChromosome()) && (
			(reg1.getStart()>= otherStart-outerRange && reg1.getStart()<= otherStop+outerRange) ||
			(reg1.getStop()>= otherStart-outerRange && reg1.getStop()<= otherStop+outerRange) ||
			(reg1.getStart()<= otherStart-outerRange && reg1.getStop()>= otherStop+outerRange) ||
			(reg1.getStart()>= otherStart-outerRange && reg1.getStop()<= otherStop+outerRange)
			);
		boolean overlap2 = chr.equals(reg2.getChromosome()) && (
			(reg2.getStart()>= otherStart-outerRange && reg2.getStart()<= otherStop+outerRange) ||
			(reg2.getStop()>= otherStart-outerRange && reg2.getStop()<= otherStop+outerRange) ||
			(reg2.getStart()<= otherStart-outerRange && reg2.getStop()>= otherStop+outerRange) ||
			(reg2.getStart()>= otherStart-outerRange && reg2.getStop()<= otherStop+outerRange)
			);
		return overlap1 || overlap2;
	}
}