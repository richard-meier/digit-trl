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

public class PositionFinder {
	private PositionRank end;
	private PositionFinder left;
	private PositionFinder right;
	private double center;
	
	public PositionFinder(double center){
		this.center = center;
		end=null;
		left=null;
		right=null;
	}

	public PositionFinder getLeft() {
		return left;
	}

	public void setLeft(PositionFinder left) {
		this.left = left;
	}

	public PositionFinder getRight() {
		return right;
	}

	public void setRight(PositionFinder right) {
		this.right = right;
	}
	
	public double getCenter(){
		return center;
	}

	public PositionRank getEnd() {
		return end;
	}

	public void setEnd(PositionRank end) {
		this.end = end;
	}
}
