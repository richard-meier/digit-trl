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

import java.util.Vector;

public class Alignment {
	public static Vector<AlignmentResult> globalNeedlemanWunsch(String a0, String b0){
		int matchWorth=2, mismatchPenalty=-2;
		int gapPenalty=-1;
		String a=a0.toUpperCase().replace("N", "n");
		String b=b0.toUpperCase();
		int[][] dpm = new int[a.length()+1][b.length()+1];
		int[][] rpm = new int[a.length()+1][b.length()+1];
		
		Vector<AlignmentStartSite> topCandidates = new Vector<AlignmentStartSite>();
		for(int i=1; i<a.length()+1; i++){
			for(int j=1; j<b.length()+1; j++){
				dpm[i][j] = dpm[i-1][j-1] + ((a.charAt(i-1)==b.charAt(j-1))?matchWorth:mismatchPenalty);
				rpm[i][j] = 0;
				
				if(dpm[i][j] < dpm[i][j-1]+gapPenalty){
					dpm[i][j] = dpm[i][j-1]+gapPenalty;
					rpm[i][j]=-1;
				}
				if(dpm[i][j] < dpm[i-1][j]+gapPenalty){
					dpm[i][j] = dpm[i-1][j]+gapPenalty;
					rpm[i][j]=1;
				}
			}
		}
		AlignmentStartSite tmp = new AlignmentStartSite();
		tmp.maxPosX=a.length();
		tmp.maxPosY=b.length();
		topCandidates.add(tmp);
		Vector<AlignmentResult> topAlignments = new Vector<AlignmentResult>();
		for(AlignmentStartSite alnStart:topCandidates){
			topAlignments.add(createAlignment(dpm,rpm,alnStart,a,b,false));
		}
		return topAlignments;
	}
	
	public static Vector<AlignmentResult> localSmithWaterman(String a0, String b0){
		int matchWorth=2, mismatchPenalty=-2;
		int gapPenalty=-1;
		String a=a0.toUpperCase().replace("N", "n");
		String b=b0.toUpperCase();
		int[][] dpm = new int[a.length()+1][b.length()+1];
		int[][] rpm = new int[a.length()+1][b.length()+1];
		
		int maxVal=Integer.MIN_VALUE;
		
		Vector<AlignmentStartSite> topCandidates = new Vector<AlignmentStartSite>();
		for(int i=1; i<a.length()+1; i++){
			for(int j=1; j<b.length()+1; j++){
				dpm[i][j] = dpm[i-1][j-1] + ((a.charAt(i-1)==b.charAt(j-1))?matchWorth:mismatchPenalty);
				rpm[i][j] = 0;
				
				if(dpm[i][j] < dpm[i][j-1]+gapPenalty){
					dpm[i][j] = dpm[i][j-1]+gapPenalty;
					rpm[i][j]=-1;
				}
				if(dpm[i][j] < dpm[i-1][j]+gapPenalty){
					dpm[i][j] = dpm[i-1][j]+gapPenalty;
					rpm[i][j]=1;
				}
				
				if(dpm[i][j]<0) dpm[i][j] = 0;
				
				if(dpm[i][j]==maxVal){
					AlignmentStartSite tmp = new AlignmentStartSite();
					tmp.maxVal=maxVal;
					tmp.maxPosX=i;
					tmp.maxPosY=j;
					topCandidates.add(tmp);
				}
				else if(dpm[i][j]>maxVal){	
					maxVal = dpm[i][j];
					topCandidates.removeAllElements();
					AlignmentStartSite tmp = new AlignmentStartSite();
					tmp.maxVal=maxVal;
					tmp.maxPosX=i;
					tmp.maxPosY=j;
					topCandidates.add(tmp);
				}
			}
		}
		Vector<AlignmentResult> topAlignments = new Vector<AlignmentResult>();
		for(AlignmentStartSite alnStart:topCandidates){
			topAlignments.add(createAlignment(dpm,rpm,alnStart,a,b,true));
		}
		return topAlignments;
	}
	
	private static void writeAlignmentScore(AlignmentResult res){
		int matchCount = 0, gapCount=0, mismatchCount=0;
		for(int i=0; i<res.alignmentString1.length(); i++){
			if(res.alignmentString1.charAt(i)==res.alignmentString2.charAt(i)){
				matchCount++;
			}
			else if(res.alignmentString1.charAt(i)=='-' || res.alignmentString2.charAt(i)=='-'){
				gapCount++;
			}
			else {
				mismatchCount++;
			}
		}
		res.score = matchCount*1.0/(1.0*(gapCount+mismatchCount+matchCount));
	}
	
	private static AlignmentResult createAlignment(int[][]dpm, int[][] rpm, AlignmentStartSite start, String a, String b, boolean smithWaterman){
		AlignmentResult output = new AlignmentResult();
		StringBuffer al1 = new StringBuffer();
		StringBuffer al2 = new StringBuffer();
		int x=start.maxPosX;
		int y=start.maxPosY;
		while(x>0 && y>0){
			if(smithWaterman && dpm[x][y]==0) break;
			if(rpm[x][y] == 0){
				al1.insert(0,a.charAt(x-1));
				al2.insert(0,b.charAt(y-1));
				x--;
				y--;
			}
			else if(rpm[x][y] == -1){
				al1.insert(0,"-");
				al2.insert(0,b.charAt(y-1));
				y--;
			}
			else if(rpm[x][y] == 1){
				al1.insert(0,a.charAt(x-1));
				al2.insert(0,"-");
				x--;
			}
		}
		output.alignmentString1=al1.toString();
		output.alignmentString2=al2.toString();
		output.startString1=a;
		output.startString2=b;
		output.sequenceStart1=x;
		output.sequenceEnd1=start.maxPosX;
		output.sequenceStart2=y;
		output.sequenceEnd2=start.maxPosY;
		writeAlignmentScore(output);
		return output;
	}
	
	public static void printMatrix(int[][] m, String a, String b){
		String a0 = "#"+a;
		String b0 = "#"+b;
		System.out.print(" \t");
		for(int j=0; j<m[0].length; j++){
			System.out.print(b0.charAt(j)+"\t");
		}
		System.out.println();
		for(int i=0; i<m.length; i++){
			System.out.print(a0.charAt(i)+"\t");
			for(int j=0; j<m[i].length; j++){
				System.out.print(m[i][j]+"\t");
			}
			System.out.println();
		}
	}
	
	private static class AlignmentStartSite{
		public int maxPosX;
		public int maxPosY;
		public int maxVal;
	}
}
