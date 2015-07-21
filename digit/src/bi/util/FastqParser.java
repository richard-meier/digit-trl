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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class FastqParser{
	private boolean hasNext;
	private String line;
	private BufferedReader reader;

	public FastqParser(File input) throws IOException{
		hasNext=true;
		reader=new BufferedReader(new FileReader(input));
		getNextLine();
	}

	public boolean hasNext(){
		return hasNext;
	}

	public FastqEntry getNext() throws IOException{
		String head=null, sequence=null, description=null, quality=null;
		do{
			if(line.charAt(0) == '@'){
				head=line;
				break;
			}
			getNextLine();
		} while(line!=null);
		getNextLine();
		sequence=line;
		getNextLine();
		if(line!=null && line.charAt(0)=='+'){
			description=line;
			getNextLine();
			quality=line;
		}
		else{
			quality=line;
		}
		getNextLine();
		if(line==null) {
			hasNext=false;
			reader.close();
		}
		return new FastqEntry(head, sequence, description, quality);
	}

	private void getNextLine() throws IOException{
		String line;
		while((line=reader.readLine())!=null){
			if(line.length()<1) continue;
			break;
		}
		this.line=line;
	}
}