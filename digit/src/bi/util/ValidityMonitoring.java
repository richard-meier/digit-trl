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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class ValidityMonitoring {
	private MappingValidation mapVal;
	private GenomeRandomFileAccess genomeAccessor;
	private static final int sampleCount=100000;
	
	public static int getSampleCount(){
		return sampleCount;
	}
	
	public ValidityMonitoring(String genomeFile, String indexFile, String chrLenFile) throws IOException{
		mapVal = new MappingValidation();
		genomeAccessor = new GenomeRandomFileAccess(genomeFile,indexFile, chrLenFile);
	}
	
	public void run(File input, File output, int outerRange) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(input));
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		String line1, line2;
		int cnt=0;
		System.out.println("GO!");
		while((line1=br.readLine())!=null){
			if(line1.length()<2) continue;
			SAMEntry read1 = new SAMEntry(line1);
			while((line2=br.readLine())!=null){
				if(line2.length()<2) continue;
				break;
			}
			if(line2 == null){
				System.err.println("ERROR: Invalid read-mate-format!");
				break;
			}
			SAMEntry read2 = new SAMEntry(line2);
			if(!read2.getQname().equals(read1.getQname())){
				System.err.println("ERROR: Ivalid format: Reads are not paired!");
				break;
			}
			String seq1 = read1.getSequence();
			String seq2 = read2.getSequence();
			if(seq1.length()<20 || seq2.length()<20) continue;
			cnt++;
			if(cnt>sampleCount) break;
			int start1 = read1.getPosition()-outerRange, start2 = read2.getPosition()-outerRange;
			int stop1 = read1.getPosition()+2*outerRange, stop2 = read2.getPosition()+2*outerRange;
			String reg1 = genomeAccessor.getGenomicSequence(read1.getRname(), start1, stop1);
			String reg2 = genomeAccessor.getGenomicSequence(read2.getRname(), start2, stop2);
			if(reg1==null || reg2==null){
				System.err.println("\tSkipping read pair!");
				continue;
			}
			
			mapVal.calculateLowerValidityScoreAndPrint(seq1, seq2, reg1, reg2, bw);
		}
		bw.close();
		br.close();
		System.out.println("DONE!");
	}
	
	public void testRun(File input, File output, int readLength) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(input));
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(output+"1P_TEST.txt"));
		String line1, line2;
		int cnt=0;
		while((line1=br.readLine())!=null){
			if(line1.length()<2) continue;
			SAMEntry read1 = new SAMEntry(line1);
			while((line2=br.readLine())!=null){
				if(line2.length()<2) continue;
				break;
			}
			if(line2 == null){
				System.err.println("ERROR: Invalid read-mate-format!");
				break;
			}
			SAMEntry read2 = new SAMEntry(line2);
			if(!read2.getQname().equals(read1.getQname())){
				System.err.println("ERROR: Ivalid format: Reads are not paired!");
				break;
			}
			String seq1 = read1.getSequence();
			String seq2 = read2.getSequence();

			int cutOff=(int)(Math.round(0.75*readLength));
			if(cutOff<50) cutOff=readLength;
			if(seq1.length()<cutOff || seq2.length()<cutOff) continue;
			
			int start1 = read1.getPosition()-readLength, start2 = read2.getPosition()-readLength;
			int stop1 = read1.getPosition()+2*readLength, stop2 = read2.getPosition()+2*readLength;
			
			if(Math.abs(read1.getPosition()-read2.getPosition())<3.01*readLength){
				continue;
			}
			cnt++;
			if(cnt%(sampleCount/10)==0){
				System.out.println("\t\t"+cnt+" pairs processed ...");
			}
			if(cnt>sampleCount) break;
			String reg1 = genomeAccessor.getGenomicSequence(read1.getRname(), start1, stop1);
			String reg2 = genomeAccessor.getGenomicSequence(read2.getRname(), start2, stop2);
			if(reg1==null || reg2==null){
				System.err.println("\tSkipping read pair!");
				continue;
			}
			
			double mvm=mapVal.calculateLowerValidityScoreAndPrint(seq1, seq2, reg1, reg2, bw);
			if(mvm<1.001 && mvm>0.999){
				bw2.write(line1);
				bw2.newLine();
				bw2.write(line2);
				bw2.newLine();
			}
		}
		bw.close();
		bw2.close();
		br.close();
	}
}
