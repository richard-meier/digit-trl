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

package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import bi.util.DnaMod;
import bi.util.FastaParser;
import bi.util.FastqEntry;
import bi.util.FastqParser;

public class ReverseComplementer {
	BufferedWriter writer;
	FastqParser fqParse;
	FastaParser faParse;
	boolean fastq;
	
	public ReverseComplementer(File input, File output, boolean fastq) throws IOException{
		this.fastq=fastq;
		if(fastq){
			fqParse = new FastqParser(input);
		}
		else{
			faParse = new FastaParser(input);
		}
		writer = new BufferedWriter(new FileWriter(output));
	}
	
	public void run() throws IOException{
		if(fastq){
			while(fqParse.hasNext()){
				FastqEntry entry = fqParse.getNext();
				writer.write(entry.getHead());
				writer.newLine();
				writer.write(DnaMod.reverseComplement(entry.getSequence()));
				writer.newLine();
				if(entry.getDescription()!=null){
					writer.write(entry.getDescription());
					writer.newLine();
				}
				writer.write(new StringBuffer(entry.getQuality()).reverse().toString());
				writer.newLine();
			}
			writer.close();
		}
		else{
			System.out.println("Using fasta instead of fastq.");
			while(faParse.hasNext()){
				FastqEntry entry = faParse.getNext();
				writer.write(entry.getHead());
				writer.newLine();
				writer.write(DnaMod.reverseComplement(entry.getSequence()));
				writer.newLine();
			}
			writer.close();
		}
		
	}
}
