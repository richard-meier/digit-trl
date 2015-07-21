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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;

public class MappingValidation {
	private HashMap<String, Double> thresholds;
	private final double NOT_MATCHING_SCORE = 0.01;
	public double lastValue=-1; 
	
	public MappingValidation(){
		thresholds = new HashMap<String, Double>();
		thresholds.put("0.1",   1.16000);
		thresholds.put("0.075", 1.18147);
		thresholds.put("0.05",  1.28395);
		thresholds.put("0.025", 1.38528);
		thresholds.put("0.01",  1.56642);
	}
	
	public MappingValidation(double targetSignificance){
		this();
		thresholds.put("custom",targetSignificance);
	}
	
	public boolean checklowerValidityScore(double lowerValidity, String pValue){
		if(!thresholds.containsKey(pValue)){
			System.err.println("ERROR::MappingValidation:checklowerValidityScore: Invalid target p-value!");
			System.exit(-1);
		}
		return thresholds.get(pValue) < lowerValidity;
	}
	
	public boolean numberOfObservationsIsLessThanTheExpectedNumberOfFalsePositives(int count, int maxCount, String pValue){
		if(!thresholds.containsKey(pValue)){
			System.err.println("ERROR::MappingValidation:checklowerValidityScore: Invalid target p-value!");
			System.exit(-1);
		}
		return numberOfObservationsIsLessThanTheExpectedNumberOfFalsePositives(count,maxCount,thresholds.get(pValue));
	}
	
	public boolean numberOfObservationsIsLessThanTheExpectedNumberOfFalsePositives(int count, int maxCount, double pValue){
		return count < maxCount*pValue;
	}
	
	public double calculateLowerValidityScore(String seq1, String seq2, String reg1, String reg2) throws IOException{
		double lv=Math.min(calculateValidityRatio(seq1,reg1,reg2), calculateValidityRatio(seq2,reg2,reg1));	
		lastValue=lv;
		return lv;
	}
	
	public double calculateLowerValidityScoreAndPrint(String seq1, String seq2, String reg1, String reg2, BufferedWriter bw) throws IOException{
		double lv = Math.min(calculateValidityRatio(seq1,reg1,reg2), calculateValidityRatio(seq2,reg2,reg1));
		if(bw !=null){
			bw.write(lv+"");
			bw.newLine();
		}
		lastValue=lv;
		return lv;
	}
	
	public double calculateLowerValidityScoreAndPrint(double validity1, double validity2, BufferedWriter bw) throws IOException{
		double lv = Math.min(validity1, validity2);
		if(bw !=null){
			bw.write(lv+"");
			bw.newLine();
		}
		lastValue=lv;
		return lv;
	}
	
	public double calculateLowerValidityScoreVERBOSE(String seq1, String seq2, String reg1, String reg2) throws IOException{
		return Math.min(calculateValidityRatioVERBOSE(seq1,reg1,reg2), calculateValidityRatioVERBOSE(seq2,reg2,reg1));	
	}
	
	private double calculateValidityRatio(String readHalf, String region1, String region2) throws IOException{
		double validityRatio=0;
		AlignmentResult tmp1 = getHighestStrandMatch(readHalf,region1);
		AlignmentResult tmp2 = getHighestStrandMatch(readHalf,region2);
		if(tmp1 != null && tmp2 != null){
			validityRatio = tmp1.lNormScore1()/tmp2.lNormScore1();
		}
		else {
			if(tmp1 != null){
				validityRatio = tmp1.lNormScore1() / NOT_MATCHING_SCORE;
			}
			else if(tmp2 != null){
				validityRatio = NOT_MATCHING_SCORE / tmp2.lNormScore1();
			}
			else{
				System.out.println("uuups: "+tmp1+"\t"+tmp2+"\t"+readHalf);
				System.out.println("WARNING:: One read half does not match at all!");
				validityRatio=1;
			}
		}
		return validityRatio;
	}
	
	private double calculateValidityRatioVERBOSE(String readHalfA, String regionA, String region2) throws IOException{
		double validityRatio=0;
		AlignmentResult tmp1 = getHighestStrandMatch(readHalfA,regionA);
		AlignmentResult tmp2 = getHighestStrandMatch(readHalfA,region2);
		if(tmp1 != null && tmp2 != null){
			validityRatio = tmp1.lNormScore1()/tmp2.lNormScore1();
		}
		else {
			if(tmp1 != null){
				validityRatio = tmp1.lNormScore1() / NOT_MATCHING_SCORE;
				System.out.println(tmp1.lNormScore1()+"\t"+NOT_MATCHING_SCORE);
				System.out.println("\t1: "+tmp1.lFrac1()+" "+tmp1.lFrac2());
			}
			else if(tmp2 != null){
				validityRatio = NOT_MATCHING_SCORE / tmp2.lNormScore1();
				System.out.println(NOT_MATCHING_SCORE+"\t"+tmp2.lNormScore1());
				System.out.println("\t2: "+tmp2.lFrac1()+" "+tmp2.lFrac2());
			}
			else{
				System.out.println("uuups: "+tmp1+"\t"+tmp2+"\t"+readHalfA);
				System.out.println("WARNING:: One read half does not match at all!");
				validityRatio=1;
			}
		}
		return validityRatio;
	}
	
	private AlignmentResult getHighestStrandMatch(String sequence, String region){
		AlignmentResult tmp=null;
		Vector<AlignmentResult> alns1 = Alignment.localSmithWaterman(sequence, region);
		Vector<AlignmentResult> alns2 = Alignment.localSmithWaterman(buildReverseComplement(sequence), region);
		double maxScore=-1;
		if(alns1.size()>0){
			for(AlignmentResult alignment:alns1){
				if(alignment.alignmentString1.length()<8) {
					continue;
				}
				double currentScore = alignment.lNormScore1();
				if(maxScore<currentScore){
					maxScore=currentScore;
					tmp=alignment;
				}
			}
		}
		if(alns2.size()>0){
			for(AlignmentResult alignment:alns2){
				if(alignment.alignmentString1.length()<8) {
					continue;
				}
				double currentScore = alignment.lNormScore1();
				if(maxScore<currentScore){
					maxScore=currentScore;
					tmp=alignment;
				}
			}
		}
		return tmp;
	}
	
	private String buildReverseComplement(String sequence){
		StringBuffer buff = new StringBuffer();
		for(int i=sequence.length()-1; i>=0; i--){
			char c = sequence.charAt(i);
			char n;
			if(c=='A'){n='T';}
			else if(c=='T'){n='A';}
			else if(c=='C'){n='G';}
			else if(c=='G'){n='C';}
			else n = 'N';
			buff.append(n);
		}
		return buff.toString();
	}
	
	public void printAndStoreReadValiditiyRatios(Vector<String>seqsMx, String region1, String region2, BufferedWriter rrw, Vector<Double> tmp) throws IOException{
		for(int ii=0; ii<seqsMx.size(); ii++){
			AlignmentResult tmp1 = getHighestStrandMatch(seqsMx.get(ii),region1);
			AlignmentResult tmp2 = getHighestStrandMatch(seqsMx.get(ii),region2);
			double frac = 0;
			if(tmp1 != null && tmp2 != null){
				frac = tmp1.lNormScore1()/tmp2.lNormScore1();
				rrw.write( (frac)+"\t");
			}
			else {
				System.out.println("uuups: "+tmp1+"\t"+tmp2+"\t"+seqsMx.get(ii));
				if(tmp1 != null){
					frac = tmp1.lNormScore1() / NOT_MATCHING_SCORE;
					rrw.write( (frac) +"\t");
				}
				else if(tmp2 != null){
					frac = NOT_MATCHING_SCORE / tmp2.lNormScore1();
					rrw.write( (frac) +"\t");
				}
			}
			tmp.add(frac);
		}
	}
}
