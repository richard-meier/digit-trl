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

public class SAMFlagCheck {
	final public static int 
		C_MULTIPLE_SEGMENTS = 			0x001,
		C_PROPERLY_ALIGNED = 			0x002,
		C_SEGMENT_UNMAPPED = 			0x004,
		C_NEXT_SEGMENT_UNMAPPED = 		0x008,
		C_REVERSED = 					0x010,
		C_NEXT_SEGMENT_REVERSED = 		0x020,
		C_FIRST_SEGEMENT = 				0x040,
		C_LAST_SEGMENT =				0x080,
		C_SECONDARY_ALIGNMENT =			0x100,
		C_NOT_PASSING_QUALITY_CONTROL = 0x200,
		C_PCR_OR_DUPLICATE = 			0x400,
		C_SUPPLEMENTARY_ALIGNMENT = 	0x800;
	
	public static boolean checkFlagGeneric(int flagValue, int flagToBeChecked){
		return (flagValue & flagToBeChecked) == flagToBeChecked;
	}
	
	public static boolean checkDublicate(int flagValue){
		return (flagValue & C_PCR_OR_DUPLICATE) == C_PCR_OR_DUPLICATE;
	}
	
	public static boolean checkReversed(int flagValue){
		return (flagValue & C_REVERSED) == C_REVERSED;
	}
	
	public static boolean checkNotPassedQualityControl(int flagValue){
		return (flagValue & C_NOT_PASSING_QUALITY_CONTROL) == C_NOT_PASSING_QUALITY_CONTROL;
	}
}
