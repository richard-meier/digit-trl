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
import java.util.Random;
import java.util.Vector;

import bi.util.Alignment;
import bi.util.AlignmentResult;

public class TestSwAlignment {

	public static void run(String outputDirectory, long testSeed) {
		System.out.println("STARTING TESTING PROCEDURE!");
		try {
			int[] testSizes = {3,5,10,20};
			System.out.println("Testing Insertions:");
			for(int i : testSizes){
				System.out.println("==> test target length of "+i);
				File outFile = null;
				if(outputDirectory != null){
					if(!outputDirectory.equals("")){
						outFile = new File(outputDirectory+".ins.sz"+i+".txt");
					}
				}
				testInsertion(testSeed, i, outFile);
			}
			System.out.println("==================================================");
			System.out.println("Testing Deletions:");
			for(int i : testSizes){
				System.out.println("==> test target length of "+i);
				File outFile = null;
				if(outputDirectory != null){
					if(!outputDirectory.equals("")){
						outFile = new File(outputDirectory+".del.sz"+i+".txt");
					}
				}
				testDeletion(testSeed, i, outFile);
			}
			System.out.println("==================================================");
			System.out.println("Testing Substitutions:");
			for(int i : testSizes){
				System.out.println("==> test target length of "+i);
				File outFile = null;
				if(outputDirectory != null){
					if(!outputDirectory.equals("")){
						outFile = new File(outputDirectory+".subs.sz"+i+".txt");
					}
				}
				testSubstitution(testSeed, i, outFile);
			}
			System.out.println("==================================================");
			System.out.println("Testing Alignment Form:");
			for(int i : testSizes){
				System.out.println("==> test target length of "+i);
				File outFile = null;
				if(outputDirectory != null){
					if(!outputDirectory.equals("")){
						outFile = new File(outputDirectory+".frm.sz"+i+".txt");
					}
				}
				testForm(testSeed, i, outFile);
			}
			
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		System.out.println("FINISHED TESTING PROCEDURE!");
	}
	
	public static void testSubstitution (long seed, int insertSize, File outFile) throws IOException{
		BufferedWriter bw = null;
		if(outFile != null){
			bw = new BufferedWriter(new FileWriter(outFile));
			bw.write("expectedSeq1\texpectedSeq2\tobservedSeq1\tobservedSeq2");
			bw.newLine();
		}
		Random myRand = new Random(seed);
		int TOTAL = 1000;
		int consensus = 0;
		int tagLength = (int)Math.round(((double)insertSize)/2); // anchor tag length
		for(int i=0; i<TOTAL; i++){
			String main   = "VVVVV"; // provide unique tail sequence
			String target = "WWWWW"; // provide unique tail sequence
			
			String front = rTemplate(myRand, insertSize);
			for(int j=0; j<tagLength; j++){front += "H";} // attach unique front tag
			
			String center1 = rTemplate(myRand, insertSize);
			String center2 = "";
			for(int j=0; j<insertSize; j++){center2 += "K";}
			
			String end = "";
			for(int j=0; j<tagLength; j++){end += "I";} // attach unique end tag
			end += rTemplate(myRand, insertSize);
			
			main += front+center1+end;
			target += front+center2+end;
			
			main   += "XXXXX"; // provide unique tail sequence
			target += "YYYYY"; // provide unique tail sequence
			
			String expected1 = front+center1+end;
			String expected2 = front+center2+end;
			
			Vector<AlignmentResult> alr = Alignment.localSmithWaterman(main, target);
			boolean found = false;
			for(int j=0; j<alr.size(); j++){
				if(
					alr.get(j).alignmentString1.equals(expected1) &&
					alr.get(j).alignmentString2.equals(expected2)
				) {
					found = true;
				}
			}
			if(found) {
				consensus++;
				if(i<3) System.out.println(alr.get(0).alignmentString1+"\n"+alr.get(0).alignmentString2);
				if(bw != null){
					
				}
			}
			if(bw != null){
				bw.write(
					expected1+"\t"+expected2+"\t"+alr.get(0).alignmentString1+"\t"+alr.get(0).alignmentString2+"\t"
				);
				bw.newLine();
			}
		}
		System.out.println("A total of "+consensus+" out of "+TOTAL+" match the expectation!");
		if(outFile != null){
			bw.close();
		}
	}
	
	public static void testDeletion (long seed, int insertSize, File outFile) throws IOException{
		BufferedWriter bw = null;
		if(outFile != null){
			bw = new BufferedWriter(new FileWriter(outFile));
			bw.write("expectedSeq1\texpectedSeq2\tobservedSeq1\tobservedSeq2");
			bw.newLine();
		}
		Random myRand = new Random(seed);
		int TOTAL = 1000;
		int consensus = 0;
		int tagLength = (int)Math.round(((double)insertSize)/2); // anchor tag length
		for(int i=0; i<TOTAL; i++){
			String main   = "VVVVV"; // provide unique tail sequence
			String target = "WWWWW"; // provide unique tail sequence
			
			String front = rTemplate(myRand,5);
			for(int j=0; j<tagLength; j++){front += "H";} // attach unique front tag
			
			String center = rTemplate(myRand,insertSize);
			
			String end = "";
			for(int j=0; j<tagLength; j++){end += "I";} // attach unique end tag
			end += rTemplate(myRand,5);
			
			main += front+center+end;
			target += front+end;
			
			main   += "XXXXX"; // provide unique tail sequence
			target += "YYYYY"; // provide unique tail sequence
			
			String expected1 = front+center+end;
			String expected2 = front;
			for(int j=0; j<insertSize; j++){ expected2 += "-";}
			expected2 += end;
			
			Vector<AlignmentResult> alr = Alignment.localSmithWaterman(main, target);
			boolean found = false;
			for(int j=0; j<alr.size(); j++){
				if(
					alr.get(j).alignmentString1.equals(expected1) &&
					alr.get(j).alignmentString2.equals(expected2)
				) {
					found = true;
				}
			}
			if(found) {
				consensus++;
				if(i<3) System.out.println(alr.get(0).alignmentString1+"\n"+alr.get(0).alignmentString2);
				if(bw != null){
					
				}
			}
			if(bw != null){
				bw.write(
					expected1+"\t"+expected2+"\t"+alr.get(0).alignmentString1+"\t"+alr.get(0).alignmentString2+"\t"
				);
				bw.newLine();
			}
		}
		System.out.println("A total of "+consensus+" out of "+TOTAL+" match the expectation!");
		if(outFile != null){
			bw.close();
		}
	}
	
	public static void testInsertion (long seed, int insertSize, File outFile) throws IOException{
		BufferedWriter bw = null;
		if(outFile != null){
			bw = new BufferedWriter(new FileWriter(outFile));
			bw.write("expectedSeq1\texpectedSeq2\tobservedSeq1\tobservedSeq2");
			bw.newLine();
		}
		Random myRand = new Random(seed);
		int TOTAL = 1000;
		int consensus = 0;
		int tagLength = (int)Math.round(((double)insertSize)/2); // anchor tag length
		for(int i=0; i<TOTAL; i++){
			String main   = "VVVVV"; // provide unique tail sequence
			String target = "WWWWW"; // provide unique tail sequence
			
			String front = rTemplate(myRand,5);
			for(int j=0; j<tagLength; j++){front += "H";} // attach unique front tag
			
			String center = rTemplate(myRand,insertSize);
			
			String end = "";
			for(int j=0; j<tagLength; j++){end += "I";} // attach unique end tag
			end += rTemplate(myRand,5);
			
			main += front+end;
			target += front+center+end;
			
			main   += "XXXXX"; // provide unique tail sequence
			target += "YYYYY"; // provide unique tail sequence
			
			String expected1 = front;
			for(int j=0; j<insertSize; j++){ expected1 += "-";}
			expected1 += end;
			String expected2 = front+center+end;
			
			Vector<AlignmentResult> alr = Alignment.localSmithWaterman(main, target);
			boolean found = false;
			for(int j=0; j<alr.size(); j++){
				if(
					alr.get(j).alignmentString1.equals(expected1) &&
					alr.get(j).alignmentString2.equals(expected2)
				) {
					found = true;
				}
			}
			if(found) {
				consensus++;
				if(i<3) System.out.println(alr.get(0).alignmentString1+"\n"+alr.get(0).alignmentString2);
			}
			if(bw != null){
				bw.write(
					expected1+"\t"+expected2+"\t"+alr.get(0).alignmentString1+"\t"+alr.get(0).alignmentString2+"\t"
				);
				bw.newLine();
			}
		}
		System.out.println("A total of "+consensus+" out of "+TOTAL+" match the expectation!");
		if(outFile != null){
			bw.close();
		}
	}
	
	public static void testForm (long seed, int insertSize, File outFile) throws IOException{
		BufferedWriter bw = null;
		if(outFile != null){
			bw = new BufferedWriter(new FileWriter(outFile));
			bw.write("valid\tobservedSeq1\tobservedSeq2");
			bw.newLine();
		}
		Random myRand = new Random(seed);
		int TOTAL = 1000;
		int consensus = 0;
		for(int i=0; i<TOTAL; i++){
			String main   = rTemplate(myRand,insertSize);
			String target = rTemplate(myRand,insertSize);
			Vector<AlignmentResult> alr = Alignment.localSmithWaterman(main, target);
			boolean found = true;
			for(int j=0; j<alr.size(); j++){
				if(alr.get(j).alignmentString1.equals("") && alr.get(j).alignmentString2.equals("")){
					continue;
				}
				boolean areCharactersValid = containsViableCharacters(alr.get(j).alignmentString1) &&
						containsViableCharacters(alr.get(j).alignmentString2);
				boolean isLengthMatching = alr.get(j).alignmentString1.length() == 
						alr.get(j).alignmentString1.length();
				boolean arePositionsViable = 
						alr.get(j).sequenceStart1 >= 0 && alr.get(j).sequenceStart1 <= main.length() &&
						alr.get(j).sequenceEnd1 >= 0 && alr.get(j).sequenceEnd1 <= main.length() &&
						alr.get(j).sequenceStart1 < alr.get(j).sequenceEnd1 &&
						alr.get(j).sequenceStart2 >= 0 && alr.get(j).sequenceStart2 <= main.length() &&
						alr.get(j).sequenceEnd2 >= 0 && alr.get(j).sequenceEnd2 <= main.length() &&
						alr.get(j).sequenceStart2 < alr.get(j).sequenceEnd2;
				if( !(areCharactersValid && isLengthMatching && arePositionsViable) ){
					found = false;
					System.out.println("ERROR: "+alr.get(j).alignmentString1);
					System.out.println("ERROR: "+alr.get(j).alignmentString2);
					break;
				}
			}
			if(found) {
				consensus++;
				if(i<3) System.out.println(alr.get(0).alignmentString1+"\n"+alr.get(0).alignmentString2);
			}
			if(bw != null){
				bw.write(
					found+"\t"+alr.get(0).alignmentString1+"\t"+alr.get(0).alignmentString2+"\t"
				);
				bw.newLine();
			}
		}
		System.out.println("A total of "+consensus+" out of "+TOTAL+" random sequences produce valid output!");
		if(outFile != null){
			bw.close();
		}
	}
	
	public static boolean containsViableCharacters(String alignmentString){
		boolean out = true;
		for(int i=0; i<alignmentString.length(); i++){
			if(
				!(
					alignmentString.charAt(i) == 'A' || alignmentString.charAt(i) == 'G' ||
					alignmentString.charAt(i) == 'C' || alignmentString.charAt(i) == 'T' ||
					alignmentString.charAt(i) == '-'
				)
			){
				out=false;
				break;
			}
		}
		return out;
	}
	
	public static char[] alph = {'A','C','T','G'};
	
	public static char rBase (Random cr) {
		return alph[cr.nextInt(4)];
	}
	
	public static String rTemplate(Random cr, int len){
		String out = "";
		for(int i=0; i<len; i++){
			out += rBase(cr);
		}
		return(out);
	}
	
}
