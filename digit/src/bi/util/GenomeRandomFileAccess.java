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
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class GenomeRandomFileAccess {
	ReentrantReadWriteLock rwl;
	private String indexFile;
	public HashMap<String, HashMap<Integer,Long>> chromosomeRoot;
	RandomAccessFile fileAccessor;
	private static final int diffStore = 10000;
	public HashMap<String,Integer> chromosomeLengthMap;
	
	public GenomeRandomFileAccess(String genomeFile, String indexFile, String chrLength) throws IOException{
		chromosomeRoot = new HashMap<String, HashMap<Integer,Long>>();
		fileAccessor = new RandomAccessFile(genomeFile,"r");
		this.indexFile=indexFile;
		if(chrLength==null){
			initialiseLengthMap();
		}
		else{
			initialiseLengthMapWithFile(new File(chrLength));
		}
		loadIndex();
		rwl = new ReentrantReadWriteLock();
	}

	public void close() throws IOException{
		fileAccessor.close();
	}
	
	public void loadIndex() throws IOException{
		File index = new File(indexFile);
		BufferedReader br = new BufferedReader(new FileReader(index));
		String line;
		String currentChromosome=null;
		while( (line=br.readLine()) != null ){
			if (line.length()==0) {
				continue;
			}
			if(line.charAt(0)=='>'){
				line=line.replace(">", "");
				currentChromosome=line;
				chromosomeRoot.put(currentChromosome, new HashMap<Integer,Long>());
				continue;
			}
			String[] entries = line.split("\\s");
			chromosomeRoot.get(currentChromosome).put(Integer.parseInt(entries[0]),Long.parseLong(entries[1]));
		}
		br.close();
	}
	
	public synchronized String getGenomicSequence(String chromosome, int targetStartPosition0, int targetStopPosition0) throws IOException{
		rwl.readLock().lock();
		try{
			int targetStopPosition=correctStop(chromosome, targetStopPosition0);
			int targetStartPosition=targetStartPosition0;
			if(targetStartPosition<1){
				targetStartPosition=1;
			}
			if(targetStartPosition>targetStopPosition){
				System.err.println("WARNING::GenomeRandomFileAccess:getGenomicSequence: Target chromosome region is out of bounds!");
				System.err.println("\t"+chromosome+":"+targetStartPosition0+"-"+targetStopPosition0+" EXCEEDS:"+chromosomeLengthMap.get(chromosome));
				return null;
			}
			StringBuffer out = new StringBuffer();
			int currentBasePosition = (targetStartPosition / diffStore)*diffStore;
			long bytePosition = chromosomeRoot.get(chromosome).get(currentBasePosition);
			fileAccessor.seek(bytePosition);
			String line;
			boolean readerIsRecordingSequence=false;
			boolean recordingIsDone=false;
			while((line=fileAccessor.readLine())!=null){
				if(line.length()==0) continue;
				if(line.charAt(0)=='>') {
					System.err.println("WARNING::GenomeRandomFileAccess:getGenomicSequence: OVERFLOW:lineskip: Target chromosome region is out of bounds!");
					System.err.println("\t"+chromosome+":"+targetStartPosition+"-"+targetStopPosition+" EXCEEDS:"+chromosomeLengthMap.get(chromosome));
					return null;
				}
				for(int i=0; i<line.length();i++){
					currentBasePosition++;
					if(currentBasePosition >= targetStopPosition){
						recordingIsDone=true;
						break;
					}
					if(currentBasePosition >= targetStartPosition){
						readerIsRecordingSequence=true;
					}
					if(readerIsRecordingSequence){
						out.append(line.charAt(i));
					}
				}
				if(recordingIsDone){
					break;
				}
			}
			if(out.toString()==null){
				System.err.println("WARNING::GenomeRandomFileAccess:getGenomicSequence: Target chromosome region is empty!");
				System.err.println("\t"+chromosome+":"+targetStartPosition+"-"+targetStopPosition);
			}
			if(out.length()-(targetStopPosition-targetStartPosition)==2){
				return out.substring(1, out.length()-1);
			}
			else if(out.length()-(targetStopPosition-targetStartPosition)==1){
				return out.substring(1, out.length());
			}
			else return out.toString();
		}
		finally{
			rwl.readLock().unlock();
		}
	}
	
	private int correctStop(String chromosome, int oldStop){
		if(oldStop > chromosomeLengthMap.get(chromosome)){
			return chromosomeLengthMap.get(chromosome);
		}
		else {
			return oldStop;
		}
	}
	
	public static void buildIndex(String genomeFile, String outFile){
		try {
			RandomAccessFile rafIn = new RandomAccessFile(genomeFile,"r");
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(outFile));
			BufferedWriter bwChr = new BufferedWriter(new FileWriter(outFile+".ChrLen.txt"));
			String line;
			int currentPosition=0;
			int lineBreakOffset=1;
			boolean first=true;
			String curChr=null;
			while((line=rafIn.readLine())!=null){
				if(line.length()==0) continue;
				if(line.charAt(0)=='>'){
					if(!first){
						bwChr.write(curChr+"\t"+currentPosition);
						bwChr.newLine();
						bwChr.flush();
					}
					first=false;
					curChr=line.substring(1);
					System.out.println("PROCESSING "+line);
					bwOut.write(line);
					bwOut.newLine();
					currentPosition=0;
					bwOut.write(currentPosition+"\t"+rafIn.getFilePointer());
					bwOut.newLine();
				}
				else{
					int addCnt=0;
					for(int i=0; i<line.length(); i++){
						addCnt++;
						currentPosition++;
						if(currentPosition%diffStore==0){
							if(i==line.length()-1){
								addCnt+=lineBreakOffset;
							}
							bwOut.write(currentPosition+"\t"+(rafIn.getFilePointer()-(line.length())+addCnt-lineBreakOffset));
							bwOut.newLine();
						}
						
						if(currentPosition%10000000==0) System.out.println("\t"+currentPosition+" reached!");
					}
				}
			}
			rafIn.close();
			bwChr.close();
			bwOut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void initialiseLengthMapWithFile(File inFile) throws IOException{
		chromosomeLengthMap = new HashMap<String,Integer>();
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		String line;
		while( (line=in.readLine())!=null ){
			if(line.length()<2) continue;
			String[] entries=line.split("\\s");
			chromosomeLengthMap.put(entries[0], Integer.parseInt(entries[1]));
		}
		in.close();
	}
	
	private void initialiseLengthMap(){
		chromosomeLengthMap = new HashMap<String,Integer>();
		chromosomeLengthMap.put("chr1", 248956422);
		chromosomeLengthMap.put("chr2", 242193529);
		chromosomeLengthMap.put("chr3", 198295559);
		chromosomeLengthMap.put("chr4", 190214555);
		chromosomeLengthMap.put("chr5", 181538259);
		chromosomeLengthMap.put("chr6", 170805979);
		chromosomeLengthMap.put("chr7", 159345973);
		chromosomeLengthMap.put("chrX", 156040895);
		chromosomeLengthMap.put("chr8", 145138636);
		chromosomeLengthMap.put("chr9", 138394717);
		chromosomeLengthMap.put("chr11", 135086622);
		chromosomeLengthMap.put("chr10", 133797422);
		chromosomeLengthMap.put("chr12", 133275309);
		chromosomeLengthMap.put("chr13", 114364328);
		chromosomeLengthMap.put("chr14", 107043718);
		chromosomeLengthMap.put("chr15", 101991189);
		chromosomeLengthMap.put("chr16", 90338345);
		chromosomeLengthMap.put("chr17", 83257441);
		chromosomeLengthMap.put("chr18", 80373285);
		chromosomeLengthMap.put("chr20", 64444167);
		chromosomeLengthMap.put("chr19", 58617616);
		chromosomeLengthMap.put("chrY", 57227415);
		chromosomeLengthMap.put("chr22", 50818468);
		chromosomeLengthMap.put("chr21", 46709983);
		chromosomeLengthMap.put("chr15_KI270905v1_alt", 5161414);
		chromosomeLengthMap.put("chr6_GL000256v2_alt", 4929269);
		chromosomeLengthMap.put("chr6_GL000254v2_alt", 4827813);
		chromosomeLengthMap.put("chr6_GL000251v2_alt", 4795265);
		chromosomeLengthMap.put("chr6_GL000253v2_alt", 4677643);
		chromosomeLengthMap.put("chr6_GL000250v2_alt", 4672374);
		chromosomeLengthMap.put("chr6_GL000255v2_alt", 4606388);
		chromosomeLengthMap.put("chr6_GL000252v2_alt", 4604811);
		chromosomeLengthMap.put("chr17_KI270857v1_alt", 2877074);
		chromosomeLengthMap.put("chr16_KI270853v1_alt", 2659700);
		chromosomeLengthMap.put("chr16_KI270728v1_random", 1872759);
		chromosomeLengthMap.put("chr17_GL000258v2_alt", 1821992);
		chromosomeLengthMap.put("chr5_GL339449v2_alt", 1612928);
		chromosomeLengthMap.put("chr14_KI270847v1_alt", 1511111);
		chromosomeLengthMap.put("chr17_KI270908v1_alt", 1423190);
		chromosomeLengthMap.put("chr14_KI270846v1_alt", 1351393);
		chromosomeLengthMap.put("chr5_KI270897v1_alt", 1144418);
		chromosomeLengthMap.put("chr7_KI270803v1_alt", 1111570);
		chromosomeLengthMap.put("chr19_GL949749v2_alt", 1091841);
		chromosomeLengthMap.put("chr19_KI270938v1_alt", 1066800);
		chromosomeLengthMap.put("chr19_GL949750v2_alt", 1066390);
		chromosomeLengthMap.put("chr19_GL949748v2_alt", 1064304);
		chromosomeLengthMap.put("chr19_GL949751v2_alt", 1002683);
		chromosomeLengthMap.put("chr19_GL949746v1_alt", 987716);
		chromosomeLengthMap.put("chr19_GL949752v1_alt", 987100);
		chromosomeLengthMap.put("chr8_KI270821v1_alt", 985506);
		chromosomeLengthMap.put("chr1_KI270763v1_alt", 911658);
		chromosomeLengthMap.put("chr6_KI270801v1_alt", 870480);
		chromosomeLengthMap.put("chr19_GL949753v2_alt", 796479);
		chromosomeLengthMap.put("chr19_GL949747v2_alt", 729520);
		chromosomeLengthMap.put("chr8_KI270822v1_alt", 624492);
		chromosomeLengthMap.put("chr4_GL000257v2_alt", 586476);
		chromosomeLengthMap.put("chr12_KI270904v1_alt", 572349);
		chromosomeLengthMap.put("chr4_KI270925v1_alt", 555799);
		chromosomeLengthMap.put("chr15_KI270852v1_alt", 478999);
		chromosomeLengthMap.put("chr15_KI270727v1_random", 448248);
		chromosomeLengthMap.put("chr9_KI270823v1_alt", 439082);
		chromosomeLengthMap.put("chr15_KI270850v1_alt", 430880);
		chromosomeLengthMap.put("chr1_KI270759v1_alt", 425601);
		chromosomeLengthMap.put("chr12_GL877876v1_alt", 408271);
		chromosomeLengthMap.put("chrUn_KI270442v1", 392061);
		chromosomeLengthMap.put("chr17_KI270862v1_alt", 391357);
		chromosomeLengthMap.put("chr15_GL383555v2_alt", 388773);
		chromosomeLengthMap.put("chr19_GL383573v1_alt", 385657);
		chromosomeLengthMap.put("chr4_KI270896v1_alt", 378547);
		chromosomeLengthMap.put("chr4_GL383528v1_alt", 376187);
		chromosomeLengthMap.put("chr17_GL383563v3_alt", 375691);
		chromosomeLengthMap.put("chr8_KI270810v1_alt", 374415);
		chromosomeLengthMap.put("chr1_GL383520v2_alt", 366580);
		chromosomeLengthMap.put("chr1_KI270762v1_alt", 354444);
		chromosomeLengthMap.put("chr15_KI270848v1_alt", 327382);
		chromosomeLengthMap.put("chr17_KI270909v1_alt", 325800);
		chromosomeLengthMap.put("chr14_KI270844v1_alt", 322166);
		chromosomeLengthMap.put("chr8_KI270900v1_alt", 318687);
		chromosomeLengthMap.put("chr10_GL383546v1_alt", 309802);
		chromosomeLengthMap.put("chr13_KI270838v1_alt", 306913);
		chromosomeLengthMap.put("chr8_KI270816v1_alt", 305841);
		chromosomeLengthMap.put("chr22_KI270879v1_alt", 304135);
		chromosomeLengthMap.put("chr8_KI270813v1_alt", 300230);
		chromosomeLengthMap.put("chr11_KI270831v1_alt", 296895);
		chromosomeLengthMap.put("chr15_GL383554v1_alt", 296527);
		chromosomeLengthMap.put("chr8_KI270811v1_alt", 292436);
		chromosomeLengthMap.put("chr18_GL383567v1_alt", 289831);
		chromosomeLengthMap.put("chrX_KI270880v1_alt", 284869);
		chromosomeLengthMap.put("chr8_KI270812v1_alt", 282736);
		chromosomeLengthMap.put("chr19_KI270921v1_alt", 282224);
		chromosomeLengthMap.put("chr17_KI270729v1_random", 280839);
		chromosomeLengthMap.put("chr17_JH159146v1_alt", 278131);
		chromosomeLengthMap.put("chrX_KI270913v1_alt", 274009);
		chromosomeLengthMap.put("chr6_KI270798v1_alt", 271782);
		chromosomeLengthMap.put("chr7_KI270808v1_alt", 271455);
		chromosomeLengthMap.put("chr22_KI270876v1_alt", 263666);
		chromosomeLengthMap.put("chr15_KI270851v1_alt", 263054);
		chromosomeLengthMap.put("chr22_KI270875v1_alt", 259914);
		chromosomeLengthMap.put("chr1_KI270766v1_alt", 256271);
		chromosomeLengthMap.put("chr19_KI270882v1_alt", 248807);
		chromosomeLengthMap.put("chr3_KI270778v1_alt", 248252);
		chromosomeLengthMap.put("chr15_KI270849v1_alt", 244917);
		chromosomeLengthMap.put("chr4_KI270786v1_alt", 244096);
		chromosomeLengthMap.put("chr12_KI270835v1_alt", 238139);
		chromosomeLengthMap.put("chr17_KI270858v1_alt", 235827);
		chromosomeLengthMap.put("chr19_KI270867v1_alt", 233762);
		chromosomeLengthMap.put("chr16_KI270855v1_alt", 232857);
		chromosomeLengthMap.put("chr8_KI270926v1_alt", 229282);
		chromosomeLengthMap.put("chr5_GL949742v1_alt", 226852);
		chromosomeLengthMap.put("chr3_KI270780v1_alt", 224108);
		chromosomeLengthMap.put("chr17_GL383565v1_alt", 223995);
		chromosomeLengthMap.put("chr2_KI270774v1_alt", 223625);
		chromosomeLengthMap.put("chr4_KI270790v1_alt", 220246);
		chromosomeLengthMap.put("chr11_KI270927v1_alt", 218612);
		chromosomeLengthMap.put("chr19_KI270932v1_alt", 215732);
		chromosomeLengthMap.put("chr11_KI270903v1_alt", 214625);
		chromosomeLengthMap.put("chr2_KI270894v1_alt", 214158);
		chromosomeLengthMap.put("chr14_GL000225v1_random", 211173);
		chromosomeLengthMap.put("chrUn_KI270743v1", 210658);
		chromosomeLengthMap.put("chr11_KI270832v1_alt", 210133);
		chromosomeLengthMap.put("chr7_KI270805v1_alt", 209988);
		chromosomeLengthMap.put("chr4_GL000008v2_random", 209709);
		chromosomeLengthMap.put("chr7_KI270809v1_alt", 209586);
		chromosomeLengthMap.put("chr19_KI270887v1_alt", 209512);
		chromosomeLengthMap.put("chr4_KI270789v1_alt", 205944);
		chromosomeLengthMap.put("chr3_KI270779v1_alt", 205312);
		chromosomeLengthMap.put("chr19_KI270914v1_alt", 205194);
		chromosomeLengthMap.put("chr19_KI270886v1_alt", 204239);
		chromosomeLengthMap.put("chr11_KI270829v1_alt", 204059);
		chromosomeLengthMap.put("chr14_GL000009v2_random", 201709);
		chromosomeLengthMap.put("chr21_GL383579v2_alt", 201197);
		chromosomeLengthMap.put("chr11_JH159136v1_alt", 200998);
		chromosomeLengthMap.put("chr19_KI270930v1_alt", 200773);
		chromosomeLengthMap.put("chrUn_KI270747v1", 198735);
		chromosomeLengthMap.put("chr18_GL383571v1_alt", 198278);
		chromosomeLengthMap.put("chr19_KI270920v1_alt", 198005);
		chromosomeLengthMap.put("chr6_KI270797v1_alt", 197536);
		chromosomeLengthMap.put("chr3_KI270935v1_alt", 197351);
		chromosomeLengthMap.put("chr17_KI270861v1_alt", 196688);
		chromosomeLengthMap.put("chr15_KI270906v1_alt", 196384);
		chromosomeLengthMap.put("chr5_KI270791v1_alt", 195710);
		chromosomeLengthMap.put("chr14_KI270722v1_random", 194050);
		chromosomeLengthMap.put("chr16_GL383556v1_alt", 192462);
		chromosomeLengthMap.put("chr13_KI270840v1_alt", 191684);
		chromosomeLengthMap.put("chr14_GL000194v1_random", 191469);
		chromosomeLengthMap.put("chr11_JH159137v1_alt", 191409);
		chromosomeLengthMap.put("chr19_KI270917v1_alt", 190932);
		chromosomeLengthMap.put("chr7_KI270899v1_alt", 190869);
		chromosomeLengthMap.put("chr19_KI270923v1_alt", 189352);
		chromosomeLengthMap.put("chr10_KI270825v1_alt", 188315);
		chromosomeLengthMap.put("chr19_GL383576v1_alt", 188024);
		chromosomeLengthMap.put("chr19_KI270922v1_alt", 187935);
		chromosomeLengthMap.put("chrUn_KI270742v1", 186739);
		chromosomeLengthMap.put("chr22_KI270878v1_alt", 186262);
		chromosomeLengthMap.put("chr19_KI270929v1_alt", 186203);
		chromosomeLengthMap.put("chr11_KI270826v1_alt", 186169);
		chromosomeLengthMap.put("chr6_KB021644v2_alt", 185823);
		chromosomeLengthMap.put("chr17_GL000205v2_random", 185591);
		chromosomeLengthMap.put("chr1_KI270765v1_alt", 185285);
		chromosomeLengthMap.put("chr19_KI270916v1_alt", 184516);
		chromosomeLengthMap.put("chr19_KI270890v1_alt", 184499);
		chromosomeLengthMap.put("chr3_KI270784v1_alt", 184404);
		chromosomeLengthMap.put("chr12_GL383551v1_alt", 184319);
		chromosomeLengthMap.put("chr20_KI270870v1_alt", 183433);
		chromosomeLengthMap.put("chrUn_GL000195v1", 182896);
		chromosomeLengthMap.put("chr1_GL383518v1_alt", 182439);
		chromosomeLengthMap.put("chr22_KI270736v1_random", 181920);
		chromosomeLengthMap.put("chr10_KI270824v1_alt", 181496);
		chromosomeLengthMap.put("chr14_KI270845v1_alt", 180703);
		chromosomeLengthMap.put("chr3_GL383526v1_alt", 180671);
		chromosomeLengthMap.put("chr13_KI270839v1_alt", 180306);
		chromosomeLengthMap.put("chr22_KI270733v1_random", 179772);
		chromosomeLengthMap.put("chrUn_GL000224v1", 179693);
		chromosomeLengthMap.put("chr10_GL383545v1_alt", 179254);
		chromosomeLengthMap.put("chrUn_GL000219v1", 179198);
		chromosomeLengthMap.put("chr5_KI270792v1_alt", 179043);
		chromosomeLengthMap.put("chr17_KI270860v1_alt", 178921);
		chromosomeLengthMap.put("chr19_GL000209v2_alt", 177381);
		chromosomeLengthMap.put("chr11_KI270830v1_alt", 177092);
		chromosomeLengthMap.put("chr9_KI270719v1_random", 176845);
		chromosomeLengthMap.put("chrUn_GL000216v2", 176608);
		chromosomeLengthMap.put("chr22_KI270928v1_alt", 176103);
		chromosomeLengthMap.put("chr1_KI270712v1_random", 176043);
		chromosomeLengthMap.put("chr6_KI270800v1_alt", 175808);
		chromosomeLengthMap.put("chr1_KI270706v1_random", 175055);
		chromosomeLengthMap.put("chr2_KI270776v1_alt", 174166);
		chromosomeLengthMap.put("chr18_KI270912v1_alt", 174061);
		chromosomeLengthMap.put("chr3_KI270777v1_alt", 173649);
		chromosomeLengthMap.put("chr5_GL383531v1_alt", 173459);
		chromosomeLengthMap.put("chr3_JH636055v2_alt", 173151);
		chromosomeLengthMap.put("chr14_KI270725v1_random", 172810);
		chromosomeLengthMap.put("chr5_KI270796v1_alt", 172708);
		chromosomeLengthMap.put("chr9_GL383541v1_alt", 171286);
		chromosomeLengthMap.put("chr19_KI270885v1_alt", 171027);
		chromosomeLengthMap.put("chr19_KI270919v1_alt", 170701);
		chromosomeLengthMap.put("chr19_KI270889v1_alt", 170698);
		chromosomeLengthMap.put("chr19_KI270891v1_alt", 170680);
		chromosomeLengthMap.put("chr19_KI270915v1_alt", 170665);
		chromosomeLengthMap.put("chr19_KI270933v1_alt", 170537);
		chromosomeLengthMap.put("chr19_KI270883v1_alt", 170399);
		chromosomeLengthMap.put("chr19_GL383575v2_alt", 170222);
		chromosomeLengthMap.put("chr19_KI270931v1_alt", 170148);
		chromosomeLengthMap.put("chr12_GL383550v2_alt", 169178);
		chromosomeLengthMap.put("chr13_KI270841v1_alt", 169134);
		chromosomeLengthMap.put("chrUn_KI270744v1", 168472);
		chromosomeLengthMap.put("chr18_KI270863v1_alt", 167999);
		chromosomeLengthMap.put("chr18_GL383569v1_alt", 167950);
		chromosomeLengthMap.put("chr12_GL877875v1_alt", 167313);
		chromosomeLengthMap.put("chr21_KI270874v1_alt", 166743);
		chromosomeLengthMap.put("chr3_KI270924v1_alt", 166540);
		chromosomeLengthMap.put("chr1_KI270761v1_alt", 165834);
		chromosomeLengthMap.put("chr3_KI270937v1_alt", 165607);
		chromosomeLengthMap.put("chr22_KI270734v1_random", 165050);
		chromosomeLengthMap.put("chr18_GL383570v1_alt", 164789);
		chromosomeLengthMap.put("chr5_KI270794v1_alt", 164558);
		chromosomeLengthMap.put("chr4_GL383527v1_alt", 164536);
		chromosomeLengthMap.put("chrUn_GL000213v1", 164239);
		chromosomeLengthMap.put("chr3_KI270936v1_alt", 164170);
		chromosomeLengthMap.put("chr3_KI270934v1_alt", 163458);
		chromosomeLengthMap.put("chr9_GL383539v1_alt", 162988);
		chromosomeLengthMap.put("chr3_KI270895v1_alt", 162896);
		chromosomeLengthMap.put("chr22_GL383582v2_alt", 162811);
		chromosomeLengthMap.put("chr3_KI270782v1_alt", 162429);
		chromosomeLengthMap.put("chr1_KI270892v1_alt", 162212);
		chromosomeLengthMap.put("chrUn_GL000220v1", 161802);
		chromosomeLengthMap.put("chr2_KI270767v1_alt", 161578);
		chromosomeLengthMap.put("chr2_KI270715v1_random", 161471);
		chromosomeLengthMap.put("chr2_KI270893v1_alt", 161218);
		chromosomeLengthMap.put("chrUn_GL000218v1", 161147);
		chromosomeLengthMap.put("chr18_GL383572v1_alt", 159547);
		chromosomeLengthMap.put("chr8_KI270817v1_alt", 158983);
		chromosomeLengthMap.put("chr4_KI270788v1_alt", 158965);
		chromosomeLengthMap.put("chrUn_KI270749v1", 158759);
		chromosomeLengthMap.put("chr7_KI270806v1_alt", 158166);
		chromosomeLengthMap.put("chr7_KI270804v1_alt", 157952);
		chromosomeLengthMap.put("chr18_KI270911v1_alt", 157710);
		chromosomeLengthMap.put("chrUn_KI270741v1", 157432);
		chromosomeLengthMap.put("chr17_KI270910v1_alt", 157099);
		chromosomeLengthMap.put("chr19_KI270884v1_alt", 157053);
		chromosomeLengthMap.put("chr19_GL383574v1_alt", 155864);
		chromosomeLengthMap.put("chr19_KI270888v1_alt", 155532);
		chromosomeLengthMap.put("chr3_GL000221v1_random", 155397);
		chromosomeLengthMap.put("chr11_GL383547v1_alt", 154407);
		chromosomeLengthMap.put("chr2_KI270716v1_random", 153799);
		chromosomeLengthMap.put("chr12_GL383553v2_alt", 152874);
		chromosomeLengthMap.put("chr6_KI270799v1_alt", 152148);
		chromosomeLengthMap.put("chr22_KI270731v1_random", 150754);
		chromosomeLengthMap.put("chrUn_KI270751v1", 150742);
		chromosomeLengthMap.put("chrUn_KI270750v1", 148850);
		chromosomeLengthMap.put("chr8_KI270818v1_alt", 145606);
		chromosomeLengthMap.put("chrX_KI270881v1_alt", 144206);
		chromosomeLengthMap.put("chr21_KI270873v1_alt", 143900);
		chromosomeLengthMap.put("chr2_GL383521v1_alt", 143390);
		chromosomeLengthMap.put("chr8_KI270814v1_alt", 141812);
		chromosomeLengthMap.put("chr12_GL383552v1_alt", 138655);
		chromosomeLengthMap.put("chrUn_KI270519v1", 138126);
		chromosomeLengthMap.put("chr2_KI270775v1_alt", 138019);
		chromosomeLengthMap.put("chr17_KI270907v1_alt", 137721);
		chromosomeLengthMap.put("chrUn_GL000214v1", 137718);
		chromosomeLengthMap.put("chr8_KI270901v1_alt", 136959);
		chromosomeLengthMap.put("chr2_KI270770v1_alt", 136240);
		chromosomeLengthMap.put("chr16_KI270854v1_alt", 134193);
		chromosomeLengthMap.put("chr8_KI270819v1_alt", 133535);
		chromosomeLengthMap.put("chr17_GL383564v2_alt", 133151);
		chromosomeLengthMap.put("chr2_KI270772v1_alt", 133041);
		chromosomeLengthMap.put("chr8_KI270815v1_alt", 132244);
		chromosomeLengthMap.put("chr5_KI270795v1_alt", 131892);
		chromosomeLengthMap.put("chr5_KI270898v1_alt", 130957);
		chromosomeLengthMap.put("chr20_GL383577v2_alt", 128386);
		chromosomeLengthMap.put("chr1_KI270708v1_random", 127682);
		chromosomeLengthMap.put("chr7_KI270807v1_alt", 126434);
		chromosomeLengthMap.put("chr5_KI270793v1_alt", 126136);
		chromosomeLengthMap.put("chr6_GL383533v1_alt", 124736);
		chromosomeLengthMap.put("chr2_GL383522v1_alt", 123821);
		chromosomeLengthMap.put("chr19_KI270918v1_alt", 123111);
		chromosomeLengthMap.put("chr12_GL383549v1_alt", 120804);
		chromosomeLengthMap.put("chr2_KI270769v1_alt", 120616);
		chromosomeLengthMap.put("chr4_KI270785v1_alt", 119912);
		chromosomeLengthMap.put("chr12_KI270834v1_alt", 119498);
		chromosomeLengthMap.put("chr7_GL383534v2_alt", 119183);
		chromosomeLengthMap.put("chr20_KI270869v1_alt", 118774);
		chromosomeLengthMap.put("chr21_GL383581v2_alt", 116689);
		chromosomeLengthMap.put("chr3_KI270781v1_alt", 113034);
		chromosomeLengthMap.put("chr17_KI270730v1_random", 112551);
		chromosomeLengthMap.put("chrUn_KI270438v1", 112505);
		chromosomeLengthMap.put("chr4_KI270787v1_alt", 111943);
		chromosomeLengthMap.put("chr18_KI270864v1_alt", 111737);
		chromosomeLengthMap.put("chr2_KI270771v1_alt", 110395);
		chromosomeLengthMap.put("chr1_GL383519v1_alt", 110268);
		chromosomeLengthMap.put("chr2_KI270768v1_alt", 110099);
		chromosomeLengthMap.put("chr1_KI270760v1_alt", 109528);
		chromosomeLengthMap.put("chr3_KI270783v1_alt", 109187);
		chromosomeLengthMap.put("chr17_KI270859v1_alt", 108763);
		chromosomeLengthMap.put("chr11_KI270902v1_alt", 106711);
		chromosomeLengthMap.put("chr18_GL383568v1_alt", 104552);
		chromosomeLengthMap.put("chr22_KI270737v1_random", 103838);
		chromosomeLengthMap.put("chr13_KI270843v1_alt", 103832);
		chromosomeLengthMap.put("chr22_KI270877v1_alt", 101331);
		chromosomeLengthMap.put("chr5_GL383530v1_alt", 101241);
		chromosomeLengthMap.put("chr11_KI270721v1_random", 100316);
		chromosomeLengthMap.put("chr22_KI270738v1_random", 99375);
		chromosomeLengthMap.put("chr22_GL383583v2_alt", 96924);
		chromosomeLengthMap.put("chr2_GL582966v2_alt", 96131);
		chromosomeLengthMap.put("chrUn_KI270748v1", 93321);
		chromosomeLengthMap.put("chrUn_KI270435v1", 92983);
		chromosomeLengthMap.put("chr5_GL000208v1_random", 92689);
		chromosomeLengthMap.put("chrUn_KI270538v1", 91309);
		chromosomeLengthMap.put("chr17_GL383566v1_alt", 90219);
		chromosomeLengthMap.put("chr16_GL383557v1_alt", 89672);
		chromosomeLengthMap.put("chr17_JH159148v1_alt", 88070);
		chromosomeLengthMap.put("chr5_GL383532v1_alt", 82728);
		chromosomeLengthMap.put("chr21_KI270872v1_alt", 82692);
		chromosomeLengthMap.put("chrUn_KI270756v1", 79590);
		chromosomeLengthMap.put("chr6_KI270758v1_alt", 76752);
		chromosomeLengthMap.put("chr12_KI270833v1_alt", 76061);
		chromosomeLengthMap.put("chr6_KI270802v1_alt", 75005);
		chromosomeLengthMap.put("chr21_GL383580v2_alt", 74653);
		chromosomeLengthMap.put("chr22_KB663609v1_alt", 74013);
		chromosomeLengthMap.put("chr22_KI270739v1_random", 73985);
		chromosomeLengthMap.put("chr9_GL383540v1_alt", 71551);
		chromosomeLengthMap.put("chrUn_KI270757v1", 71251);
		chromosomeLengthMap.put("chr2_KI270773v1_alt", 70887);
		chromosomeLengthMap.put("chr17_JH159147v1_alt", 70345);
		chromosomeLengthMap.put("chr11_KI270827v1_alt", 67707);
		chromosomeLengthMap.put("chr1_KI270709v1_random", 66860);
		chromosomeLengthMap.put("chrUn_KI270746v1", 66486);
		chromosomeLengthMap.put("chr16_KI270856v1_alt", 63982);
		chromosomeLengthMap.put("chr21_GL383578v2_alt", 63917);
		chromosomeLengthMap.put("chrUn_KI270753v1", 62944);
		chromosomeLengthMap.put("chr19_KI270868v1_alt", 61734);
		chromosomeLengthMap.put("chr9_GL383542v1_alt", 60032);
		chromosomeLengthMap.put("chr20_KI270871v1_alt", 58661);
		chromosomeLengthMap.put("chr12_KI270836v1_alt", 56134);
		chromosomeLengthMap.put("chr19_KI270865v1_alt", 52969);
		chromosomeLengthMap.put("chr1_KI270764v1_alt", 50258);
		chromosomeLengthMap.put("chrUn_KI270589v1", 44474);
		chromosomeLengthMap.put("chr14_KI270726v1_random", 43739);
		chromosomeLengthMap.put("chr19_KI270866v1_alt", 43156);
		chromosomeLengthMap.put("chr22_KI270735v1_random", 42811);
		chromosomeLengthMap.put("chr1_KI270711v1_random", 42210);
		chromosomeLengthMap.put("chrUn_KI270745v1", 41891);
		chromosomeLengthMap.put("chr1_KI270714v1_random", 41717);
		chromosomeLengthMap.put("chr22_KI270732v1_random", 41543);
		chromosomeLengthMap.put("chr1_KI270713v1_random", 40745);
		chromosomeLengthMap.put("chrUn_KI270754v1", 40191);
		chromosomeLengthMap.put("chr1_KI270710v1_random", 40176);
		chromosomeLengthMap.put("chr12_KI270837v1_alt", 40090);
		chromosomeLengthMap.put("chr9_KI270717v1_random", 40062);
		chromosomeLengthMap.put("chr14_KI270724v1_random", 39555);
		chromosomeLengthMap.put("chr9_KI270720v1_random", 39050);
		chromosomeLengthMap.put("chr14_KI270723v1_random", 38115);
		chromosomeLengthMap.put("chr9_KI270718v1_random", 38054);
		chromosomeLengthMap.put("chrUn_KI270317v1", 37690);
		chromosomeLengthMap.put("chr13_KI270842v1_alt", 37287);
		chromosomeLengthMap.put("chrY_KI270740v1_random", 37240);
		chromosomeLengthMap.put("chrUn_KI270755v1", 36723);
		chromosomeLengthMap.put("chr8_KI270820v1_alt", 36640);
		chromosomeLengthMap.put("chr1_KI270707v1_random", 32032);
		chromosomeLengthMap.put("chrUn_KI270579v1", 31033);
		chromosomeLengthMap.put("chrUn_KI270752v1", 27745);
		chromosomeLengthMap.put("chrUn_KI270512v1", 22689);
		chromosomeLengthMap.put("chrUn_KI270322v1", 21476);
		chromosomeLengthMap.put("chrM", 16569);
		chromosomeLengthMap.put("chrUn_GL000226v1", 15008);
		chromosomeLengthMap.put("chrUn_KI270311v1", 12399);
		chromosomeLengthMap.put("chrUn_KI270366v1", 8320);
		chromosomeLengthMap.put("chrUn_KI270511v1", 8127);
		chromosomeLengthMap.put("chrUn_KI270448v1", 7992);
		chromosomeLengthMap.put("chrUn_KI270521v1", 7642);
		chromosomeLengthMap.put("chrUn_KI270581v1", 7046);
		chromosomeLengthMap.put("chrUn_KI270582v1", 6504);
		chromosomeLengthMap.put("chrUn_KI270515v1", 6361);
		chromosomeLengthMap.put("chrUn_KI270588v1", 6158);
		chromosomeLengthMap.put("chrUn_KI270591v1", 5796);
		chromosomeLengthMap.put("chrUn_KI270522v1", 5674);
		chromosomeLengthMap.put("chrUn_KI270507v1", 5353);
		chromosomeLengthMap.put("chrUn_KI270590v1", 4685);
		chromosomeLengthMap.put("chrUn_KI270584v1", 4513);
		chromosomeLengthMap.put("chrUn_KI270320v1", 4416);
		chromosomeLengthMap.put("chrUn_KI270382v1", 4215);
		chromosomeLengthMap.put("chrUn_KI270468v1", 4055);
		chromosomeLengthMap.put("chrUn_KI270467v1", 3920);
		chromosomeLengthMap.put("chrUn_KI270362v1", 3530);
		chromosomeLengthMap.put("chrUn_KI270517v1", 3253);
		chromosomeLengthMap.put("chrUn_KI270593v1", 3041);
		chromosomeLengthMap.put("chrUn_KI270528v1", 2983);
		chromosomeLengthMap.put("chrUn_KI270587v1", 2969);
		chromosomeLengthMap.put("chrUn_KI270364v1", 2855);
		chromosomeLengthMap.put("chrUn_KI270371v1", 2805);
		chromosomeLengthMap.put("chrUn_KI270333v1", 2699);
		chromosomeLengthMap.put("chrUn_KI270374v1", 2656);
		chromosomeLengthMap.put("chrUn_KI270411v1", 2646);
		chromosomeLengthMap.put("chrUn_KI270414v1", 2489);
		chromosomeLengthMap.put("chrUn_KI270510v1", 2415);
		chromosomeLengthMap.put("chrUn_KI270390v1", 2387);
		chromosomeLengthMap.put("chrUn_KI270375v1", 2378);
		chromosomeLengthMap.put("chrUn_KI270420v1", 2321);
		chromosomeLengthMap.put("chrUn_KI270509v1", 2318);
		chromosomeLengthMap.put("chrUn_KI270315v1", 2276);
		chromosomeLengthMap.put("chrUn_KI270302v1", 2274);
		chromosomeLengthMap.put("chrUn_KI270518v1", 2186);
		chromosomeLengthMap.put("chrUn_KI270530v1", 2168);
		chromosomeLengthMap.put("chrUn_KI270304v1", 2165);
		chromosomeLengthMap.put("chrUn_KI270418v1", 2145);
		chromosomeLengthMap.put("chrUn_KI270424v1", 2140);
		chromosomeLengthMap.put("chrUn_KI270417v1", 2043);
		chromosomeLengthMap.put("chrUn_KI270508v1", 1951);
		chromosomeLengthMap.put("chrUn_KI270303v1", 1942);
		chromosomeLengthMap.put("chrUn_KI270381v1", 1930);
		chromosomeLengthMap.put("chrUn_KI270529v1", 1899);
		chromosomeLengthMap.put("chrUn_KI270425v1", 1884);
		chromosomeLengthMap.put("chrUn_KI270396v1", 1880);
		chromosomeLengthMap.put("chrUn_KI270363v1", 1803);
		chromosomeLengthMap.put("chrUn_KI270386v1", 1788);
		chromosomeLengthMap.put("chrUn_KI270465v1", 1774);
		chromosomeLengthMap.put("chrUn_KI270383v1", 1750);
		chromosomeLengthMap.put("chrUn_KI270384v1", 1658);
		chromosomeLengthMap.put("chrUn_KI270330v1", 1652);
		chromosomeLengthMap.put("chrUn_KI270372v1", 1650);
		chromosomeLengthMap.put("chrUn_KI270548v1", 1599);
		chromosomeLengthMap.put("chrUn_KI270580v1", 1553);
		chromosomeLengthMap.put("chrUn_KI270387v1", 1537);
		chromosomeLengthMap.put("chrUn_KI270391v1", 1484);
		chromosomeLengthMap.put("chrUn_KI270305v1", 1472);
		chromosomeLengthMap.put("chrUn_KI270373v1", 1451);
		chromosomeLengthMap.put("chrUn_KI270422v1", 1445);
		chromosomeLengthMap.put("chrUn_KI270316v1", 1444);
		chromosomeLengthMap.put("chrUn_KI270338v1", 1428);
		chromosomeLengthMap.put("chrUn_KI270340v1", 1428);
		chromosomeLengthMap.put("chrUn_KI270583v1", 1400);
		chromosomeLengthMap.put("chrUn_KI270334v1", 1368);
		chromosomeLengthMap.put("chrUn_KI270429v1", 1361);
		chromosomeLengthMap.put("chrUn_KI270393v1", 1308);
		chromosomeLengthMap.put("chrUn_KI270516v1", 1300);
		chromosomeLengthMap.put("chrUn_KI270389v1", 1298);
		chromosomeLengthMap.put("chrUn_KI270466v1", 1233);
		chromosomeLengthMap.put("chrUn_KI270388v1", 1216);
		chromosomeLengthMap.put("chrUn_KI270544v1", 1202);
		chromosomeLengthMap.put("chrUn_KI270310v1", 1201);
		chromosomeLengthMap.put("chrUn_KI270412v1", 1179);
		chromosomeLengthMap.put("chrUn_KI270395v1", 1143);
		chromosomeLengthMap.put("chrUn_KI270376v1", 1136);
		chromosomeLengthMap.put("chrUn_KI270337v1", 1121);
		chromosomeLengthMap.put("chrUn_KI270335v1", 1048);
		chromosomeLengthMap.put("chrUn_KI270378v1", 1048);
		chromosomeLengthMap.put("chrUn_KI270379v1", 1045);
		chromosomeLengthMap.put("chrUn_KI270329v1", 1040);
		chromosomeLengthMap.put("chrUn_KI270419v1", 1029);
		chromosomeLengthMap.put("chrUn_KI270336v1", 1026);
		chromosomeLengthMap.put("chrUn_KI270312v1", 998);
		chromosomeLengthMap.put("chrUn_KI270539v1", 993);
		chromosomeLengthMap.put("chrUn_KI270385v1", 990);
		chromosomeLengthMap.put("chrUn_KI270423v1", 981);
		chromosomeLengthMap.put("chrUn_KI270392v1", 971);
		chromosomeLengthMap.put("chrUn_KI270394v1", 970);
	}
}
