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

public interface PositionIdentifier {
	
	public PositionIdentifier getPrevious();
	public void setPrevious(PositionIdentifier previous);
	public PositionIdentifier getNext();
	public void setNext(PositionIdentifier next);
	public int getPosition();
	public ReadTag getTag();
	public String getChromosomeName();
	public PositionIdentifier getMate();
	public void setMate(PositionIdentifier mate);
}
