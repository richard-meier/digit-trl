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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bi.util.AnnotationCheck;
import bi.util.SAMEntry;
import bi.util.SAMFlagCheck;

public class Analyser {
	private String outputPath;
	private File inputFile;
	private int readCounter=0;
	private double standardDeviation=0;
	private double medianValueOfReadLength;
	private int insCount, delCount, transCount;
	private double numberTimesSigma;
	private boolean inputIsBWA, inputIsBowtie;
	private Pattern patternSoftclippedMatchSoftclipped = Pattern.compile("^\\d+[S](\\d+[MDI]+)+\\d+[S]$");
	private Pattern patternMatchSoftclipped = Pattern.compile("^(\\d+[MDI]{1})+(\\d+[S]){1}$");
	private Pattern patternSoftclippedMatch = Pattern.compile("^(\\d+[S]){1}(\\d+[MDI]{1})+$");
	private Pattern patternSoftclipped = Pattern.compile("\\d+[S]");
	private Pattern patternCntOfSoftclipped = Pattern.compile("\\d+");
	private AnnotationCheck lowComplexityAnnotationCheck,testAC;

	private int inversionCounter;
	private double mapQual=0;

	public Analyser(String inputFileStart, String outputPath, double numberTimesSigma, String type, File lca_file, double mapQual, String chrLen) throws IOException{
		this.mapQual=mapQual;
		this.inputFile=new File(inputFileStart);
		this.outputPath=outputPath;
		this.numberTimesSigma=numberTimesSigma;
		inversionCounter=0;
		if(type !=null && type.contains("bowtie")){
			inputIsBowtie=true;
			System.out.println("Assuming alignment has been created with bowtie2");
		}
		else {
			inputIsBWA=true;
			System.out.println("Assuming alignment has been created with bwa");
		}
		lowComplexityAnnotationCheck = new AnnotationCheck(lca_file,chrLen);
		System.out.println("Detected MAPQ-TH of "+this.mapQual);
	}

	public void run() throws IOException{
		System.out.println("Doing weird things... -__- 0");
		System.out.println("Calculating Statistical Parameters");
		calulateMedian();
		System.out.println("\tStep 1/2 done...");
		estimateStandardDeviation();
		System.out.println("\tStep 2/2 done.");
		System.out.println("Applying Filters");
		produceInDelTrans();
		System.out.println("Done");
	}

	private void calulateMedian() throws IOException{
		BufferedReader brs = new BufferedReader(new FileReader(inputFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputPath+"MAPQ-Dist.txt")));
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(new File(outputPath+"SEP-Dist.txt")));
		String sCurrentLineStart;
		String sCurrentLineEnd;
		Vector<Integer> allLengthValues = new Vector<Integer>();
		int cnt=0;
		int switchCount=0;
		while (((sCurrentLineStart = brs.readLine()) != null)) {
			if(sCurrentLineStart.length()<2) continue;
			if(sCurrentLineStart.charAt(0)=='@') continue;
			sCurrentLineEnd = brs.readLine();
			if(sCurrentLineEnd == null){
				System.out.println("WARNING:: The last read of the input file is not paired!");
				break;
			}
			cnt++;
			SAMEntry r1 = new SAMEntry(sCurrentLineStart);
			SAMEntry r2 = new SAMEntry(sCurrentLineEnd);
			
			if(cnt<1000000){
				bw.write(r1.getMappingQuality()+""); bw.newLine();
				bw.write(r2.getMappingQuality()+""); bw.newLine();
			}
			
			byte file1=SAMFlagCheck.checkFlagGeneric(r1.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			byte file2=SAMFlagCheck.checkFlagGeneric(r2.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			boolean swtch=false;
			if(file1==(byte)2 && file2==(byte)1){
				SAMEntry tmp = r1;
				r1=r2;
				r2=tmp;
				swtch=true;
			}
			if(!r1.getQname().equals(r2.getQname())){
				System.err.println("ERROR:: There exist subsequent unpaired reads in the input file!\nSort the sam file according to read names and re-run the program!");
				System.exit(-1);
			}
			if(cnt<1000000){
				bw2.write((r1.getPosition()-r2.getPosition())+""); bw2.newLine();
			}
			if(areEntriesIrrelevantBWA(r1, r2)) continue;
			if(swtch) switchCount++;
			int length=Math.abs(r1.getPosition()-r2.getPosition());
			allLengthValues.add(length);
		}
		System.out.println("\t"+allLengthValues.size()+" of "+cnt+" pairs used!");
		System.out.println("\t"+switchCount+" readOrigins had to be switched!");
		Collections.sort(allLengthValues);
		readCounter=allLengthValues.size();
		if(readCounter%2==1)
			medianValueOfReadLength=allLengthValues.get((readCounter+1)/2);
		else{
			if(readCounter!=0){
				int value1=allLengthValues.get(readCounter/2);
				int value2=allLengthValues.get((readCounter/2)+1);
				medianValueOfReadLength=(value1+value2)/2.0;
			}else{
				System.out.println("EMPTY LIST !!!");
				medianValueOfReadLength=0;
			}
		}
		System.out.println("\tcalculated median size = "+medianValueOfReadLength);
		brs.close();
		bw.close();
		bw2.close();
	}

	private void estimateStandardDeviation() throws IOException{
		BufferedReader brs = new BufferedReader(new FileReader(inputFile));
		String sCurrentLineStart;
		String sCurrentLineEnd;
		int currentReadLength;
		double differenceToMean;
		int windowCounter=0;
		while (((sCurrentLineStart = brs.readLine()) != null)) {
			if(sCurrentLineStart.length()<2) continue;
			if(sCurrentLineStart.charAt(0)=='@') continue;
			sCurrentLineEnd = brs.readLine();
			if(sCurrentLineEnd == null){
				System.out.println("WARNING:: The last read of the input file is not paired!");
				break;
			}
			SAMEntry r1 = new SAMEntry(sCurrentLineStart);
			SAMEntry r2 = new SAMEntry(sCurrentLineEnd);
			byte file1=SAMFlagCheck.checkFlagGeneric(r1.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			byte file2=SAMFlagCheck.checkFlagGeneric(r2.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			if(file1==(byte)2 && file2==(byte)1){
				SAMEntry tmp = r1;
				r1=r2;
				r2=tmp;
			}
			if(!r1.getQname().equals(r2.getQname())){
				System.err.println("ERROR:: There exist subsequent unpaired reads in the input file!\nSort the sam file according to read names and re-run the program!");
				System.exit(-1);
			}
			if(areEntriesIrrelevantBWA(r1, r2)) continue;
			if(Math.abs(r1.getPosition()-r2.getPosition()) <= medianValueOfReadLength*2) windowCounter++;
		}
		brs.close();
		brs = new BufferedReader(new FileReader(inputFile));
		while (((sCurrentLineStart = brs.readLine()) != null)) {
			if(sCurrentLineStart.length()<2) continue;
			if(sCurrentLineStart.charAt(0)=='@') continue;
			sCurrentLineEnd = brs.readLine();
			if(sCurrentLineEnd == null){
				System.out.println("WARNING:: The last read of the input file is not paired!");
				break;
			}
			SAMEntry r1 = new SAMEntry(sCurrentLineStart);
			SAMEntry r2 = new SAMEntry(sCurrentLineEnd);
			if(!r1.getQname().equals(r2.getQname())){
				System.err.println("ERROR:: There exist subsequent unpaired reads in the input file!\nSort the sam file according to read names and re-run the program!");
				System.exit(-1);
			}
			if(areEntriesIrrelevantBWA(r1, r2)) continue;
			if(Math.abs(r1.getPosition()-r2.getPosition()) <= medianValueOfReadLength*2){
				currentReadLength=Math.abs(r1.getPosition()-r2.getPosition());
				differenceToMean=currentReadLength-medianValueOfReadLength;
				standardDeviation+=differenceToMean*differenceToMean/(windowCounter-1.0);
			}
		}
		brs.close();
		standardDeviation=Math.sqrt(standardDeviation);
	}

	private void produceInDelTrans() throws IOException{
		BufferedReader brs = new BufferedReader(new FileReader(inputFile));
		String currentFileName = removeDotSAM(inputFile.getName());
		if(currentFileName==null) currentFileName = outputPath;
		else currentFileName = outputPath+currentFileName+"_";
		BufferedWriter bwin = new BufferedWriter(new FileWriter(currentFileName+"ins.sam"));
		BufferedWriter bwd = new BufferedWriter(new FileWriter(currentFileName+"del.sam"));
		BufferedWriter bwt = new BufferedWriter(new FileWriter(currentFileName+"trl.sam"));
		BufferedWriter bwiv = new BufferedWriter(new FileWriter(currentFileName+"inv.sam"));
		BufferedWriter bwunmR1 = new BufferedWriter(new FileWriter(currentFileName+"unmapped_R1.sam"));
		BufferedWriter bwunmR2 = new BufferedWriter(new FileWriter(currentFileName+"unmapped_R2.sam"));
		BufferedWriter badQualR1 = new BufferedWriter(new FileWriter(currentFileName+"badQ_R1.sam"));
		BufferedWriter badQualR2 = new BufferedWriter(new FileWriter(currentFileName+"badQ_R2.sam"));
		BufferedWriter concordantWriter = new BufferedWriter(new FileWriter(currentFileName+"conc.sam"));
		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(currentFileName+"summary.txt"));
		BufferedWriter softsummary = new BufferedWriter(new FileWriter(currentFileName+"softsummary.fq"));
		BufferedWriter softclips = new BufferedWriter(new FileWriter(currentFileName+"softclips.sam"));
		BufferedWriter lowComp = new BufferedWriter(new FileWriter(currentFileName+"low_complexity.sam"));
		
		BufferedWriter bwDM = new BufferedWriter(new FileWriter(outputPath+"DISC_MAPQ.txt"));
		String line1;
		String line2;
		int currentReadLength;
		int XOR_DUP_COUNT = 0;
		int AND_DUP_COUNT = 0;
		int TOO_SHORT_COUNT = 0;
		while (((line1 = brs.readLine()) != null)) {
			if(line1.length()<2) continue;
			if(line1.charAt(0)=='@') continue;
			line2 = brs.readLine();
			if(line2 == null){
				System.out.println("WARNING:: The last read of the input file is not paired!");
				break;
			}
			SAMEntry read1 = new SAMEntry(line1);
			SAMEntry read2 = new SAMEntry(line2);
			byte file1=SAMFlagCheck.checkFlagGeneric(read1.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			byte file2=SAMFlagCheck.checkFlagGeneric(read2.getFlag(), SAMFlagCheck.C_FIRST_SEGEMENT) ? (byte)1 :(byte)2;
			if(file1==(byte)2 && file2==(byte)1){
				SAMEntry tmp = read1;
				read1=read2;
				read2=tmp;
			}
			if(!read1.getQname().equals(read2.getQname())){
				System.err.println("ERROR:: There exist subsequent unpaired reads in the input file!\nSort the sam file according to read names and re-run the program!");
				System.exit(-1);
			}
			
			if(entryDoesNotMap(read1)|| entryDoesNotMap(read2)) {
				continue;
			}
			String splitter = ":;|;:";
			if(matchAndSoftClipOnly(read1,read2)){
				boolean found = true;
				if(matchAndSoftClipOnly(read1) ){
					if(getClipedReadSeqAndQuality(read1)!=null){ 
						String temp = getClipedReadSeqAndQuality(read1);
						String[] tmp = temp.split("z{1}!{1}~{1}#{1}z{1}!{1}z{1}#{1}~{1}!{1}z{1}"); 
						softsummary.write("@"+read1.getQname()+splitter+read1.getCigar()+splitter);
						softsummary.write(read2.getRname()+"_"+read2.getPosition()+splitter);
						softsummary.write(read1.getRname()+"_"+read1.getPosition()+splitter);
						softsummary.write(tmp[1]);
						softsummary.newLine();
						softsummary.write(tmp[0]);
						softsummary.newLine();
					}
				}else if(matchAndSoftClipOnly(read2) ){
					if(getClipedReadSeqAndQuality(read2)!=null){ 
						String temp = getClipedReadSeqAndQuality(read2);
						String[] tmp = temp.split("z{1}!{1}~{1}#{1}z{1}!{1}z{1}#{1}~{1}!{1}z{1}"); 
						softsummary.write("@"+read1.getQname()+splitter+read2.getCigar()+splitter);
						softsummary.write(read1.getRname()+"_"+read1.getPosition()+splitter);
						softsummary.write(read2.getRname()+"_"+read2.getPosition()+splitter);
						softsummary.write(tmp[1]);
						softsummary.newLine();
						softsummary.write(tmp[0]);
						softsummary.newLine();
					}
				} else {
					found=false;
				}
				if(found){
					softclips.write(line1);
					softclips.newLine();
					softclips.write(line2);
					softclips.newLine();
				}
				continue;
			}
			if(SAMFlagCheck.checkDublicate(read1.getFlag()) ^ SAMFlagCheck.checkDublicate(read2.getFlag())){
				XOR_DUP_COUNT++;
			}
			if(
					( hasBadQuality(read1) || hasBadQuality(read2) ) ||
					( Math.abs(read1.getPosition()-read2.getPosition())<5 && read1.getRname().equals(read2.getRname()) ) ||
					( SAMFlagCheck.checkDublicate(read1.getFlag()) && SAMFlagCheck.checkDublicate(read2.getFlag()) ) ||
					( read1.getSequence().length() < 20 || read2.getSequence().length() < 20)
			) {
				if(SAMFlagCheck.checkDublicate(read1.getFlag()) && SAMFlagCheck.checkDublicate(read2.getFlag())){
					AND_DUP_COUNT++;
				}
				if(SAMFlagCheck.checkDublicate(read1.getFlag()) && SAMFlagCheck.checkDublicate(read2.getFlag())){
					TOO_SHORT_COUNT++;
				}
				continue;
			}

			if(
					checkLowComplexity(read1) || checkLowComplexity(read2)
			){
				continue;
			}
			
			if (entriesRepresentTranslocation(read1, read2)){
				transCount++;
				bwt.write(line1+"\n"+line2+"\n\n");
			} else{
				if(entriesRepresentInversion(read1, read2)){
					bwDM.write(read1.getMappingQuality()+""); bwDM.newLine();
					bwiv.write(line1+"\n"+line2+"\n\n");
					inversionCounter++;
				}
				else{
					currentReadLength=Math.abs(read1.getPosition()-read2.getPosition());
					if(currentReadLength<(medianValueOfReadLength-(numberTimesSigma*standardDeviation))){ //sigma=standardDeviation
						insCount++;
						bwin.write(line1+"\n"+line2+"\n\n");
					}
					else if (currentReadLength>(medianValueOfReadLength+(numberTimesSigma*standardDeviation))){
						delCount++;
						bwd.write(line1+"\n"+line2+"\n\n");
					}
					else{
						if(testAC != null){
							if(testAC.getOverlap(read1.getRname(), read1.getPosition(), read1.getPosition()+100)>0){
								continue;
							}
						}
						
						if(
							!( read1.getMappingQuality()<10 || read2.getMappingQuality()<10 )
						){
							concordantWriter.write(line1);
							concordantWriter.newLine();
							concordantWriter.write(line2);
							concordantWriter.newLine();
						}
						
					}
				}
			}
		}
		System.out.println(XOR_DUP_COUNT+" reads have been wrongly assigned by picard tools as single strand dublications!");
		System.out.println(AND_DUP_COUNT+" reads have been filtered out because they are marked as paired dublications!");
		System.out.println(TOO_SHORT_COUNT+" reads have been filtered out because they are too short!");
		String summary = "";
		summary += "Insertions: "+insCount+"\n"+"Deletions: "+delCount+"\n"+"Translocations: "+transCount+"\n"+"Potential Inversions: "+inversionCounter;
		summary += "\nstandardDeviation: "+standardDeviation;
		summary += "\nmedian: "+medianValueOfReadLength;
		summary += "\nnumber of read-pairs: "+readCounter;
		int threshold = (int)(Math.floor( 2.33 * standardDeviation + medianValueOfReadLength)) + 1;
		summary += "\nrecommended cluster threshold: "+threshold;
		summaryWriter.write(summary);
		System.out.println(summary);
		bwin.close();
		bwd.close();
		bwt.close();
		brs.close();
		bwunmR1.close();
		bwunmR2.close();
		badQualR1.close();
		badQualR2.close();
		bwiv.close();
		concordantWriter.close();
		summaryWriter.close();
		softsummary.close();
		softclips.close();
		lowComp.close();
		bwDM.close();
	}

	final static String myPat = "z!~#z!z#~!z";
	private String getClipedReadSeqAndQuality(SAMEntry read) {
		String cigar = read.getCigar();
		Matcher matcherMatchSoftclipped = patternSoftclippedMatchSoftclipped.matcher(cigar);
		if(matcherMatchSoftclipped.find())return null;
		
		matcherMatchSoftclipped = patternMatchSoftclipped.matcher(cigar);
		Matcher matcherSoftclippedMatch = patternSoftclippedMatch.matcher(cigar);
		
		if(matcherMatchSoftclipped.find()){
			Matcher matcherSoftclipped = patternSoftclipped.matcher(matcherMatchSoftclipped.group());
			if(matcherSoftclipped.find()){
				Matcher matcherCntOfSoftclipped = patternCntOfSoftclipped.matcher(matcherSoftclipped.group());
				if(matcherCntOfSoftclipped.find()){
					int softClipedCnt = Integer.parseInt(matcherCntOfSoftclipped.group());
					if(softClipedCnt<20) return null;
					String returnString = read.getSequence().substring(read.getSequence().length()-softClipedCnt);
					returnString+="\n+\n";
					returnString+=read.getQuality().substring(read.getSequence().length()-softClipedCnt);
					returnString+=myPat+read.getSequence().substring(0,read.getSequence().length()-softClipedCnt);
					return returnString;
				}else{
					System.err.println("SoftClipPatter: Should not happen! Missing number?");
					System.exit(-1);
				}
			}else{
				System.err.println("SoftClipPatter: Should not happen!");
				System.exit(-1);
			}
		}else if(matcherSoftclippedMatch.find()){
			Matcher matcherSoftclipped = patternSoftclipped.matcher(matcherSoftclippedMatch.group());
			if(matcherSoftclipped.find()){
				Matcher matcherCntOfSoftclipped = patternCntOfSoftclipped.matcher(matcherSoftclipped.group());
				if(matcherCntOfSoftclipped.find()){
					int softClipedCnt = Integer.parseInt(matcherCntOfSoftclipped.group());
					if(softClipedCnt<20) return null;
					String returnString = read.getSequence().substring(0,softClipedCnt);
					returnString+="\n+\n";
					returnString+=read.getQuality().substring(0,softClipedCnt);
					returnString+=myPat+read.getSequence().substring(softClipedCnt,read.getSequence().length());
					return returnString;
				}else{
					System.err.println("SoftClipPatter: Should not happen! Missing number?");
					System.exit(-1);
				}
			}else{
				System.err.println("SoftClipPatter: Should not happen!");
				System.exit(-1);
			}
		}else {
			return null;
		}
		return null;
	}

	private boolean checkLowComplexity(SAMEntry read){
		boolean doesOverlap = lowComplexityAnnotationCheck.getOverlap(read.getRname(), read.getPosition(), read.getPosition()+100)>9;
		if(doesOverlap){
				return true;
		}
		else 
			return false;
	}
	
	private boolean matchAndSoftClipOnly(SAMEntry read1, SAMEntry read2) {
		if(matchAndSoftClipOnly(read1) || matchAndSoftClipOnly(read2)) return true;
		return false;
	}

	private boolean matchAndSoftClipOnly(SAMEntry read) {
		Pattern pattern = Pattern.compile("\\d+[S]+");
		Matcher matcher = pattern.matcher(read.getCigar());
		if(matcher.find()) return true;
		return false;
	}

	private String removeDotSAM(String line){
		if(line.length()<5) return null;
		if(line.charAt(line.length()-1)=='m' &&
				line.charAt(line.length()-2)=='a' &&
				line.charAt(line.length()-3)=='s' &&
				line.charAt(line.length()-4)=='.'
				) return line.substring(0, line.length()-4);
		else return line;
	}

	private boolean areEntriesIrrelevantBWA(SAMEntry r1, SAMEntry r2) {
		if(hasBadQuality(r1) || hasBadQuality(r2)) {
			return true;
		}
		if(Math.abs(r1.getPosition()-r2.getPosition())<5 && r1.getRname().equals(r2.getRname())) {
			return true;
		}
		if(entriesRepresentTranslocation(r1,r2)) {
			return true;
		}
		if(entryDoesNotMap(r1) || entryDoesNotMap(r2)) {
			return true;
		}
		if(entriesRepresentInversion(r1,r2)) {
			return true;
		}
		return false;
	}

	private boolean entriesRepresentTranslocation(SAMEntry r1, SAMEntry r2){
		String chromosomeName1 = r1.getRname();
		String chromosomeName2 = r2.getRname();
		return !(chromosomeName1.equals(chromosomeName2));
	}

	private boolean entryDoesNotMap(SAMEntry re){
		return re.getPosition() == 0;
	}

	private boolean entriesRepresentInversion(SAMEntry r1, SAMEntry r2){
		int flagR1 = r1.getFlag();
		int flagR2 = r2.getFlag();
		if(entryStrandsAlternate(r1,r2)){
			if(!SAMFlagCheck.checkReversed(flagR1) && SAMFlagCheck.checkReversed(flagR2)){ // inversion check
				if(r2.getPosition()-r1.getPosition()<0) return true;
			}
			else if(SAMFlagCheck.checkReversed(flagR1) && !SAMFlagCheck.checkReversed(flagR2)){
				if(r1.getPosition()-r2.getPosition()<0) return true;
			}
		}
		return false;
	}

	private boolean entryStrandsAlternate(SAMEntry r1, SAMEntry r2){
		return SAMFlagCheck.checkReversed(r1.getFlag()) ^ SAMFlagCheck.checkReversed(r2.getFlag());
	}

	private boolean hasBadQuality(SAMEntry e){
		if(e.getSequence().length()<20) return true;
		if(e.getMappingQuality()<mapQual || SAMFlagCheck.checkNotPassedQualityControl(e.getFlag())) return true;
		else if(inputIsBWA){
			if(!e.getXT().equals("XT:A:U")) return true;
			else return false;
		}
		else return false;
	}

	public String[] getRevChrPosBWA(String line){
		String[] Entries=line.split("\\s");
		String[] revChrPos={Entries[1],Entries[2], Entries[3]};
		return revChrPos;
	}
}

