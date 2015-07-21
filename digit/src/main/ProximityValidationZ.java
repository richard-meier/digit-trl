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
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bi.util.AllPurposeNode;
import bi.util.Cluster;
import bi.util.GenomeRandomFileAccess;
import bi.util.MappingValidation;
import bi.util.MateCluster;
import bi.util.ProximityCheckZ;
import bi.util.ReadTag;
import bi.util.SAMEntry;
import bi.util.SAMFlagCheck;
import bi.util.StringPair;
import bi.util.ValidityMonitoring;
import bi.util.ZIndependentMvmFilter;

public class ProximityValidationZ {
	private MappingValidation mapVal;
	private static final double P_VALUE = 0.001;
	
	private GenomeRandomFileAccess genomeAccessor;
	
	private Vector<Integer> positionTrackerForTemporalCalculations;
	private Vector<String> chromosomeNames;
	private ProximityCheckZ prox;
	private File inFile, outFile;
	private Vector<ReadTag> emptyList1, emptyList2;
	private int proximityCount,proximityThreshold;
	private static final long[] magnitudes = { 
		1000l, 10000l, 100000l, 1000000l, 10000000l, 100000000l, 1000000000l, 10000000000l, 100000000000l, 1000000000000l, 
		10000000000000l, 100000000000000l, 1000000000000000l, 10000000000000000l, 100000000000000000l,1000000000000000000l
	};
	
	private String genomeFile, indexFile, chrLenFile;
 	
	private int numberOfReadsInFile;
	private double pairedValidityRatioThreshold;
	
	public ProximityValidationZ(File input, File output, int proximityCount, int proximityThreshold, String genomeFile, String indexFile, boolean buildTree, String chrLenFile) throws IOException{
		this.chrLenFile=chrLenFile;
		positionTrackerForTemporalCalculations=new Vector<Integer>();
		genomeAccessor = new GenomeRandomFileAccess(genomeFile, indexFile, chrLenFile);
		numberOfReadsInFile=0;
		this.inFile=input;
		this.outFile=output;
		this.genomeFile=genomeFile;
		this.indexFile=indexFile;
		
		prox = new ProximityCheckZ(proximityThreshold,Integer.MAX_VALUE/10);
		this.proximityCount = proximityCount;
		this.proximityThreshold = proximityThreshold;
		
		chromosomeNames = new Vector<String>();
		emptyList1 = new Vector<ReadTag>();
		emptyList2 = new Vector<ReadTag>();
		
		if(buildTree){
			System.out.println("Initialising Search Super Structure");
			recordRefNames(inFile);
			initialiseSearchTree(inFile,false);
		}
	}
	
	private void reinitialise(boolean onlyTakePerfectMatches) throws IOException{
		numberOfReadsInFile=0;
		prox = new ProximityCheckZ(proximityThreshold,Integer.MAX_VALUE/10);
		chromosomeNames = new Vector<String>();
		emptyList1 = new Vector<ReadTag>();
		emptyList2 = new Vector<ReadTag>();
		System.out.println("Initialising Search Super Structure");
		recordRefNames(inFile);
		initialiseSearchTree(inFile,onlyTakePerfectMatches);
	}
	
	private void recordRefNames(File file) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line="";
		int magnitudeSurpassIndex = 0;
		System.out.println("\tScanning File For Chromosomes");
		long startingTime = System.currentTimeMillis();
		while((line=br.readLine())!=null){
			if(line.length()<2) continue;
			SAMEntry r1 = new SAMEntry(line);
			while((line=br.readLine())!=null){
				if(line.length()<2) continue;
				break;
			}
			if(line == null){
				break;
			}
			SAMEntry r2 = new SAMEntry(line);
			numberOfReadsInFile++;
			if(numberOfReadsInFile > magnitudes[magnitudeSurpassIndex]){
				System.out.println("\t\t"+magnitudes[magnitudeSurpassIndex]+" reads processed [time = "+(System.currentTimeMillis()-startingTime)+"ms] ...");
				magnitudeSurpassIndex++;
			}
			if(!chromosomeNames.contains(transformChromosomeName(r1.getRname()))) chromosomeNames.add(transformChromosomeName(r1.getRname()));
			if(!chromosomeNames.contains(transformChromosomeName(r2.getRname()))) chromosomeNames.add(transformChromosomeName(r2.getRname()));
		}
		System.out.println("\t\tTotal of "+numberOfReadsInFile+" reads in file [time = "+(System.currentTimeMillis()-startingTime)+"ms]");
		br.close();
	}
	
	private void initialiseSearchTree(File file, boolean onlyTakePerfectMatches) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		String line2;
		int magnitudeSurpassIndex = 0;
		System.out.println("\tBuilding Search Links");
		long startingTime = System.currentTimeMillis();
		int readCount=0;
		while((line=br.readLine())!=null){
			if(line.length()<2) continue;
			SAMEntry r1 = new SAMEntry(line);
			while((line2=br.readLine())!=null){
				if(line2.length()<2) continue;
				break;
			}
			if(line2 == null){
				break;
			}
			SAMEntry r2 = new SAMEntry(line2);
			if(onlyTakePerfectMatches){
				Pattern pattern1 = Pattern.compile("\\d+[ISD]+");
				Matcher matcher1 = pattern1.matcher(r1.getCigar());
				if (matcher1.find()){
					continue;
				}
				Pattern pattern2 = Pattern.compile("\\d+[ISD]+");
				Matcher matcher2 = pattern2.matcher(r2.getCigar());
				if (matcher2.find()){
					continue;
				}
			}
			readCount ++;
			if(numberOfReadsInFile > 1000){
				if(readCount % (numberOfReadsInFile/20) == 0){
					double perc = Math.round(100.0*readCount/numberOfReadsInFile);
					System.out.println("\t\t"+perc+"% processed [time = "+(System.currentTimeMillis()-startingTime)+"ms] ...");
				}
			}
			if(readCount > magnitudes[magnitudeSurpassIndex]){
				System.out.println("\t\t"+magnitudes[magnitudeSurpassIndex]+" reads processed [time = "+(System.currentTimeMillis()-startingTime)+"ms] ...");
				magnitudeSurpassIndex++;
			}
			char strand1=SAMFlagCheck.checkReversed(r1.getFlag()) ? '-' : '+';
			char strand2=SAMFlagCheck.checkReversed(r2.getFlag()) ? '-' : '+';
			prox.add(new ReadTag(r1.getQname(),(byte)1), new ReadTag(r2.getQname(),(byte)2), r1.getPosition(), r2.getPosition(), 
				transformChromosomeName(r1.getRname()), transformChromosomeName(r2.getRname()), strand1, strand2);
		}
		br.close();
		System.out.println("\tTransformimg tree leafs into linked list ...");
		prox.finishTree(readCount);
		System.out.println("\t\t100% processed [time = "+(System.currentTimeMillis()-startingTime)+"ms]");
	}
	
	public Vector<Cluster> run(boolean printReadNames, boolean onlyTakePerfectMatches) throws IOException{
		File folder = new File(outFile.getParent());
		File[] listOfFiles = folder.listFiles();
		File summaryFile=null;
		for(File f:listOfFiles){
			if(f.getAbsolutePath().contains("_summary.txt")){
				summaryFile = f;
				break;
			}
		}
		double numberOfReadPairs = 15215665;
		
		BufferedReader br_summary = new BufferedReader(new FileReader(summaryFile));
		String lin;
		lin=br_summary.readLine(); lin=br_summary.readLine(); lin=br_summary.readLine(); lin=br_summary.readLine(); lin=br_summary.readLine();
		lin=br_summary.readLine(); lin=br_summary.readLine();
		br_summary.close();
		numberOfReadPairs = Double.parseDouble(lin.split("\\s")[3]);
		System.out.println(numberOfReadPairs);
		RD_CNT=0;
		int test=0;
		int maxClusterSize=-1;
		System.out.println("Building Proximity Clusters");
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		System.out.println("out2 = "+outFile.getAbsolutePath()+".circos.txt");
		BufferedWriter bwcirc = new BufferedWriter(new FileWriter(new File(outFile.getAbsolutePath()+".circos.txt")));
		HashMap<String,String> readTagNameToLines = new HashMap<String,String>();

		String line1;
		String line2;
		Vector<ReadTag> alreadyUsed = new Vector<ReadTag>();
		Vector<Cluster> clusterPairs = new Vector<Cluster>();
		Vector<Integer> sizes = new Vector<Integer>();
		while((line1=br.readLine())!=null){
			if(line1.length()<2) continue;
			SAMEntry r1 = new SAMEntry(line1);
			while((line2=br.readLine())!=null){
				if(line2.length()<2) continue;
				break;
			}
			if(line2 == null){
				System.err.println("ERROR: Invalid read-mate-format");
				break;
			}
			SAMEntry r2 = new SAMEntry(line2);
			if(onlyTakePerfectMatches){
				Pattern pattern1 = Pattern.compile("\\d+[ISD]+");
				Matcher matcher1 = pattern1.matcher(r1.getCigar());
				if (matcher1.find()){
					continue;
				}
				Pattern pattern2 = Pattern.compile("\\d+[ISD]+");
				Matcher matcher2 = pattern2.matcher(r2.getCigar());
				if (matcher2.find()){
					continue;
				}
			}
			emptyList1.removeAllElements();
			emptyList2.removeAllElements();
			boolean currentReadIsAlreadyPartOfCluster=alreadyUsed.contains(prox.getTagByName(r1.getQname(), 1));
			ReadTag currentTag1 = prox.getTagByName(r1.getQname(),1);
			ReadTag currentTag2 = prox.getTagByName(r1.getQname(),2);
			
			readTagNameToLines.put(currentTag1.getName()+currentTag1.getFileNumber(), line1);
			readTagNameToLines.put(currentTag2.getName()+currentTag2.getFileNumber(), line2);
			boolean clusterGotFiltered=false;
			if(!currentReadIsAlreadyPartOfCluster){
				prox.findClusterNeighbours(emptyList1, emptyList2, r1.getQname());
				int start1=Integer.MAX_VALUE; int stop1=Integer.MIN_VALUE;
				int start2=Integer.MAX_VALUE; int stop2=Integer.MIN_VALUE;
				for(ReadTag t:emptyList1){
					int pos=prox.getAllPurposeNode(t.getName()+t.getFileNumber()).getPosition();
					if(start1>pos) start1=pos;
					if(stop1<pos) stop1=pos;
				}
				for(ReadTag t:emptyList2){
					int pos=prox.getAllPurposeNode(t.getName()+t.getFileNumber()).getPosition();
					if(start2>pos) start2=pos;
					if(stop2<pos) stop2=pos;
				}
			}
			if(emptyList1.size()>=proximityCount || currentReadIsAlreadyPartOfCluster){
				RD_CNT+=1;
				if(!currentReadIsAlreadyPartOfCluster && !clusterGotFiltered){
					alreadyUsed.addAll(emptyList1);
					alreadyUsed.addAll(emptyList2);
					sizes.add(emptyList1.size());
					addNewClusters(clusterPairs, emptyList1,emptyList2, bwcirc);
				}else if(!clusterGotFiltered){
					
				}
				if(maxClusterSize<emptyList1.size()) maxClusterSize = emptyList1.size();
				emptyList1.removeAllElements();
				emptyList2.removeAllElements();
			}
			test++;
			if(numberOfReadsInFile>100 && test%(numberOfReadsInFile/10)==0) System.out.println(" --> "+(Math.round(100.0*test/numberOfReadsInFile))+"% ("+test+") max_clr="+maxClusterSize);
		}
		br.close();
		Collections.sort(sizes);
		int readCounter=sizes.size();
		double medianValueOfReadLength=0;
		if(readCounter==1){
			medianValueOfReadLength=sizes.get(0);
		}
		else if(readCounter==2){
			medianValueOfReadLength=(sizes.get(0)+sizes.get(1))/2.0;
		}
		else if(readCounter%2==1)
			medianValueOfReadLength=sizes.get((readCounter+1)/2);
		else{
			if(readCounter!=0){
				int value1=sizes.get(readCounter/2);
				int value2=sizes.get((readCounter/2)+1);
				medianValueOfReadLength=(value1+value2)/2.0;
			}else{
				medianValueOfReadLength=0;
			}
		}
		double sum=0;
		for(int i:sizes){
			sum+=i;
		}
		double mean = sum/(1.0*sizes.size());
		double sum2=0;
		for(int i:sizes){
			sum2+=Math.abs(i-mean);
		}
		double stdev = sum2/(1.0*sizes.size());
		System.out.println("median cluster count: "+medianValueOfReadLength);
		System.out.println("mean cluster count: "+mean);
		System.out.println("stdev cluster count: "+stdev);

		for(int i=0; i<clusterPairs.size();i++){
			Cluster cl = clusterPairs.get(i);
			MateCluster cluster1 = cl.getMate1();
			MateCluster cluster2 = cl.getMate2();
			StringBuffer circosOutput = new StringBuffer();
			if(cluster1.getChromosomeName().equals("hsM") || cluster2.getChromosomeName().equals("hsM")) continue;;
			cluster1.appendCircosEntryTo(circosOutput);
			cluster2.appendCircosEntryTo(circosOutput);
			bwcirc.write(circosOutput.toString());
			bwcirc.newLine();
			circosOutput.delete(0, circosOutput.length());
			circosOutput.append("#\t"+cluster1.getOccupants().size()+"\t");
			cluster1.appendOccupantNames(circosOutput);
			bwcirc.write(circosOutput.toString());
			bwcirc.newLine();
			circosOutput.delete(0, circosOutput.length());
			circosOutput.append("#\t"+cluster2.getOccupants().size()+"\t");
			cluster2.appendOccupantNames(circosOutput);
			bwcirc.write(circosOutput.toString());
			bwcirc.newLine();
			
			for(ReadTag rt:cl.getMate1().getOccupants()){
				ReadTag rt2=null;
				for(ReadTag rtt:cl.getMate2().getOccupants()){
					if(rtt.getName().equals(rt.getName())){
						rt2=rtt;
						break;
					}
				}
				
				bw.write(readTagNameToLines.get(rt.getName()+rt.getFileNumber()));
				bw.newLine();
				bw.write(readTagNameToLines.get(rt2.getName()+rt2.getFileNumber()));
				bw.newLine();
			}
		}

		bw.close();
		bwcirc.close();
		System.out.println("Done");
		return clusterPairs;
	}

	private void filterReadPairsAccordingToMvm(BufferedWriter logW) throws IOException{
		System.out.println("Start filtering");
		BufferedWriter test = new BufferedWriter(new FileWriter(new File(outFile.getParent()+"/filtered_cluster_reads.txt")));
		BufferedWriter pairs = new BufferedWriter(new FileWriter(new File(outFile.getParent()+"/pairs_remapped_c"+proximityCount+".txt")));
		BufferedWriter rejected = new BufferedWriter(new FileWriter(new File(outFile.getParent()+"/rejected_pairs_c"+proximityCount+".txt")));
		BufferedWriter approved = new BufferedWriter(new FileWriter(new File(outFile.getParent()+"/approved_pairs_c"+proximityCount+".txt")));
		
		BufferedReader br = new BufferedReader(new FileReader(outFile.getAbsolutePath()));
		String line1, line2;
		int allCount=0, filterCount=0, accessErrorCount=0;;
		while((line1=br.readLine())!=null){
			if(line1.length()<2) continue;
			SAMEntry r1 = new SAMEntry(line1);
			while((line2=br.readLine())!=null){
				if(line2.length()<2) continue;
				break;
			}
			if(line2 == null){
				System.err.println("ERROR: Invalid read-mate-format");
				break;
			}
			SAMEntry r2 = new SAMEntry(line2);
			allCount++;
			if(RD_CNT>100 && allCount%(RD_CNT/10)==0){
				System.out.println(" --> "+(Math.round(100.0*allCount/RD_CNT))+"%");
			}
			
			if(doesValidityFilterFail(r1,r2,proximityThreshold,pairs)){
				if(mapVal.lastValue<0) {
					accessErrorCount++;
					continue;
				}
				filterCount++;
				rejected.write(r1.getRname()+":"+r1.getPosition()+"\t"+r2.getRname()+":"+r2.getPosition()+"\t"+mapVal.lastValue);
				rejected.newLine();
				continue;
			}
			else{
				approved.write(r1.getRname()+":"+r1.getPosition()+"\t"+r2.getRname()+":"+r2.getPosition()+"\t"+mapVal.lastValue);
				approved.newLine();
			}
			if(line1.contains("\tZM:Z:")){
				test.write(line1);
				test.newLine();
				test.write(line2);
				test.newLine();
			}
			else{
				test.write(line1+"\tZM:Z:"+mapVal.lastValue);
				test.newLine();
				test.write(line2+"\tZM:Z:"+mapVal.lastValue);
				test.newLine();
			}
		}
		System.out.println("Filtering done: \n\t\t"+filterCount+"/"+allCount+" reads removed");
		logW.write("\t"+filterCount+"/"+allCount+" reads removed");
		logW.newLine();
		System.out.println("\t\t"+accessErrorCount+"/"+allCount+" reads skipped due to invalid mapping positions (access errors)");
		pairs.close();
		test.close();
		rejected.close();
		approved.close();
		br.close();
	}
	
	private void independentMvmFilterMultiThread(File in, File out, int threadNumber) throws IOException, InterruptedException{
		System.out.println("MVM filter: Start reading input pairs");
		ConcurrentLinkedQueue<StringPair> input = new ConcurrentLinkedQueue<StringPair>();
		BufferedReader br = new BufferedReader(new FileReader(in));
		String line;
		long totalProgress=0;
		while( (line=br.readLine())!=null ){
			if(line.length()<2) continue;
			if(line.charAt(0)=='@') continue;
			String line2;
			while( (line2=br.readLine())!=null ){
				if(line2.length()<2) continue;
				if(line2.charAt(0)=='@') continue;
				break;
			}
			
			SAMEntry e1 = new SAMEntry(line);
			SAMEntry e2 = new SAMEntry(line2);
			
			if(!e1.getQname().equals(e2.getQname())){
				System.err.println("ERROR:ProximityValidation:independentMvmFilterMultiThread: Unpaired reads in input file!");
				System.err.println("\t"+e1.getQname());
				System.err.println("\t"+e2.getQname());
				System.exit(-1);
			}
			totalProgress++;
			
			input.add(new StringPair(line,line2));
		}
		br.close();
		System.out.println("\tTotal of "+totalProgress+" reads processed!");
		
		System.out.println("MVM filter: Start filtering using "+threadNumber+" threads");
		ZIndependentMvmFilter[] pool = new ZIndependentMvmFilter[threadNumber];
		Thread[] threads = new Thread[threadNumber];
		File[] tempFiles = new File[threadNumber];
		String folder = !out.isDirectory() ? out.getParent() : out.getAbsolutePath();
		if(folder.charAt(folder.length()-1)!='/') folder += '/';
		for(int t=0; t<threadNumber; t++){
			tempFiles[t]=new File(folder+"thread"+t+"_temp.txt");
			pool[t]=new ZIndependentMvmFilter(input, tempFiles[t], genomeAccessor, pairedValidityRatioThreshold, proximityThreshold, RLEN);
			threads[t]=new Thread(pool[t]);
			threads[t].start();
			System.out.println("\tThread"+t+" running!");
		}
		
		int parts=50;
		int amount = 100/parts;
		long percGap = totalProgress/parts;
		int percLim=1;
		boolean jobIsStillRunning=true;
		
		while(jobIsStillRunning){
			int finishedThreads=0;
			jobIsStillRunning=false;
			int sum=0;
			for(int t=0; t<threadNumber; t++){
				sum+=pool[t].getProgress();
				if(threads[t].isAlive()) {
					jobIsStillRunning=true;
				}
				else{
					finishedThreads++;
				}
			}
			if(totalProgress>1000 && sum>=percLim*percGap){
				System.out.println(" --> "+(percLim*amount)+"% --> "+finishedThreads+" threads finished.");
				percLim+=1;
			}
			Thread.sleep(1000);
		}
		
		for(int t=0; t<threadNumber; t++){
			threads[t].join();
		}
		
		System.out.println("MVM filter: Start combining thread temp files");
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		for(int t=0; t<threadNumber; t++){
			BufferedReader rd = new BufferedReader(new FileReader(tempFiles[t]));
			while( (line=rd.readLine())!=null ){
				bw.write(line);
				bw.newLine();
			}
			rd.close();
		}
		bw.close();
		for(int t=0; t<threadNumber; t++){
			if(tempFiles[t].exists()){
				tempFiles[t].delete();
			}
		}
		
		System.out.println("MVM filter: FINISHED!");
	}
	
	private static int RD_CNT=0;
	private int RLEN=-1;
	
	public void performMvmFilter(File in, File out, int readLength) throws IOException, InterruptedException{
		System.out.println("Calculating read validity threshold");
		pairedValidityRatioThreshold=obtainPairedValidityThreshold(P_VALUE,readLength);
		System.out.println("\tThreshold="+pairedValidityRatioThreshold);
		RLEN=readLength;
		mapVal = new MappingValidation(pairedValidityRatioThreshold);
		independentMvmFilterMultiThread(in, out, 7);
	}
	
	public void performClusterAnalysis(int readLength) throws IOException{
		BufferedWriter logW = new BufferedWriter(new FileWriter(new File(outFile.getParent()+"/prox_log.txt")));
		System.out.println("Calculating read validity threshold");
		pairedValidityRatioThreshold=obtainPairedValidityThreshold(P_VALUE,readLength);
		System.out.println("\tThreshold="+pairedValidityRatioThreshold);
		logW.write("threshold="+pairedValidityRatioThreshold);
		logW.newLine();
		mapVal = new MappingValidation(pairedValidityRatioThreshold);
		RLEN=readLength;
		Vector<Cluster> translocations;
		File in1 = inFile;
		File out1 = outFile;
		System.out.println("out = "+outFile.getAbsolutePath());
		translocations = run(true,false);
		System.out.println("\tNumber of clusters before filtering: "+translocations.size());
		logW.write("number of clusters before filtering: "+translocations.size());
		logW.newLine();
		filterReadPairsAccordingToMvm(logW);
		inFile = new File(in1.getParent()+"/filtered_cluster_reads.txt");
		outFile = new File(out1.getAbsolutePath()+".reclustered.txt");
		System.out.println("out = "+outFile.getAbsolutePath());
		
		reinitialise(false);
		translocations = run(true,false);
		System.out.println("\tNumber of clusters after filtering: "+translocations.size());
		logW.write("number of clusters after filtering: "+translocations.size());
		logW.newLine();
		inFile = outFile;
		outFile = new File(out1.getAbsolutePath()+".brPts.txt");
		int remainingCount=0;
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out1.getAbsolutePath()+".circos.txt")));
		for(Cluster translocation:translocations){
			StringBuffer out = new StringBuffer();
			if(doesFoundClusterConsistOfEnoughUniqueEvents(translocation)){
				appendTranslocationOutputFormat(out, translocation, "-");
				bw.write(out.toString());
				bw.newLine();
				remainingCount++;
			}
		}
		System.out.println("\tNumber of clusters after pcr uniqueness filter: "+remainingCount);
		logW.write("number of clusters after pcr uniqueness filter: "+remainingCount);
		logW.newLine();
		bw.close();
		logW.close();
	}
	
	private boolean doesFoundClusterConsistOfEnoughUniqueEvents(Cluster translocation){
		Vector<ReadTag> reads1 = translocation.getMate1().getOccupants();
		Vector<ReadTag> reads2 = translocation.getMate2().getOccupants();
		boolean eventsUnique = checkUniqueEvents(reads1) && checkUniqueEvents(reads2);
		boolean sizeIsTrustworthy = translocation.getLengthSum()*0.5 >= proximityCount;
		return eventsUnique && sizeIsTrustworthy;
	}
	
	private boolean checkUniqueEvents(Vector<ReadTag> reads){
		positionTrackerForTemporalCalculations.removeAllElements();
		for(ReadTag read:reads){
			int pos = read.getPositionTag().getPosition();
			if(!positionTrackerForTemporalCalculations.contains(pos)){
				positionTrackerForTemporalCalculations.add(pos);
			}
		}
		return positionTrackerForTemporalCalculations.size() >= proximityCount;
	}
	
	private double obtainPairedValidityThreshold(double significanceLevel, int readLength) throws IOException{
		File concordantFile=grepConcFile();
		if(concordantFile==null){
			System.err.println("ERROR::ProximityValidation:obtainPairedValidityThreshold: Specified input file is in a folder that does not contain a concordant-read-file!");
			System.exit(-1);
		}
		File concordantRatioFile = new File(concordantFile.getParent()+"/pairs_remapped_conc.txt");
		ValidityMonitoring valmo = new ValidityMonitoring(genomeFile,indexFile, chrLenFile);
		valmo.run(concordantFile, concordantRatioFile, readLength);
		Vector<Double> values = new Vector<Double>();
		BufferedReader br = new BufferedReader(new FileReader(concordantRatioFile));
		String line1;
		while( (line1=br.readLine()) != null ){
			if(line1.length()<1) continue;
			values.add(Double.parseDouble(line1));
		}
		br.close();
		Collections.sort(values);
		double threshold = getThresholdFromCdf(values,0.005);
		return threshold;
	}
	
	private File grepConcFile(){
		File folder = inFile.getParentFile();
		File[] listOfFiles = folder.listFiles();
		for(File f:listOfFiles){
			if(f.getName().contains("conc.sam")){
				return f;
			}
		}
		return null;
	}
	
	private static double getThresholdFromCdf(Vector<Double>data, double propThreshold){
		double out=data.get(0);
		for(int i=1; i<data.size();i++){
			double prop = getSimplifiedEdfForSortedValues(data.size(), i);
			if(prop>propThreshold){
				out=data.get(i-1);
				break;
			}
		}
		return out;
	}
	
	private static double getSimplifiedEdfForSortedValues(double n, double j){
		if(j>=n) return 1;
		else return j/n;
	}
		
	private void appendTranslocationOutputFormat(StringBuffer out, Cluster translocation, String type){
		translocation.getMate1().appendCircosEntryTo(out);
		translocation.getMate2().appendCircosEntryTo(out);
		out.append("\n");
		out.append("#\t"+translocation.getMate1().getOccupants().size()+"\t");
		translocation.getMate1().appendOccupantNames(out);
		out.append("\n");
		out.append("#\t"+translocation.getMate2().getOccupants().size()+"\t");
		translocation.getMate2().appendOccupantNames(out);
		out.append("\n");
	}
	
	private void addNewClusters(Vector<Cluster> pairs, Vector<ReadTag> emptyList1, Vector<ReadTag> emptyList2, BufferedWriter bwc) throws IOException {
		MateCluster cluster1 = getNewCluster(emptyList1);
		MateCluster cluster2 = getNewCluster(emptyList2);
		cluster1.setMate(cluster2); cluster1.addOccupants(emptyList1);
		cluster2.setMate(cluster1); cluster2.addOccupants(emptyList2);
		pairs.add(new Cluster(cluster1,cluster2));		
	}

	private MateCluster getNewCluster(Vector<ReadTag> emptyList) {
		int start=Integer.MAX_VALUE;
		int end=Integer.MIN_VALUE;
		String chrName = prox.getAllPurposeNodeByName(emptyList.get(0).getName(), emptyList.get(0).getFileNumber()).getChromosomeName();
		for(ReadTag tag : emptyList){
			AllPurposeNode rank = prox.getAllPurposeNodeByName(tag.getName(),tag.getFileNumber());
			if(rank.getPosition()<start) start = rank.getPosition();
			if(rank.getPosition()>end) end = rank.getPosition();
		}
		return new MateCluster(start,end,chrName);
	}
	
	private static String transformChromosomeName(String name){
		if(name.charAt(0)=='c' && name.charAt(1)=='h'&& name.charAt(2)=='r'){
			return name.replace("chr", "hs");
		}
		else return name;
	}
	
	private boolean doesValidityFilterFail(SAMEntry read1, SAMEntry read2, int outerRange, BufferedWriter bw) throws IOException{
		if(read1==null || read2==null){
			System.out.println("READS ARE NULL");
			System.out.println(read1);
			System.out.println(read2);
		}
		
		String seq1 = read1.getSequence();
		String seq2 = read2.getSequence();
		int start1 = read1.getPosition()-outerRange, start2 = read2.getPosition()-outerRange;
		int stop1 = read1.getPosition()+outerRange, stop2 = read2.getPosition()+outerRange;
		boolean regionCoordinatesOverlap = 
				(start1<=start2 && stop1>=start2) ||
				(start1<=stop2 && stop1>=stop2) ||
				(start1>=start2 && stop1 <= stop2);
		boolean extendedClustersOverlap = read1.getRname().equals(read2.getRname()) && regionCoordinatesOverlap;
		if(extendedClustersOverlap){
			if(Math.abs(start1-stop2)>Math.abs(start2-stop1)){
				int diff=Math.abs(read2.getPosition()-(read1.getPosition()+RLEN))/2 - 1;
				stop1=read1.getPosition()+RLEN+diff;
				start2=read2.getPosition()-diff;
			} else{
				int diff=Math.abs(read1.getPosition()-(read2.getPosition()+RLEN))/2 - 1;
				stop2=read2.getPosition()+RLEN+diff;
				start1=read1.getPosition()-diff;
			}
		}
		String reg1 = genomeAccessor.getGenomicSequence(read1.getRname(), start1, stop1);
		String reg2 = genomeAccessor.getGenomicSequence(read2.getRname(), start2, stop2);
		double lowerValidity;
		if(reg1==null || reg2==null){
			mapVal.lastValue=-1;
			return true;
		}
		
		lowerValidity = mapVal.calculateLowerValidityScoreAndPrint(seq1, seq2, reg1, reg2, bw);
		if(mapVal.checklowerValidityScore(lowerValidity, "custom")) return false;
		else return true;
	}
}
