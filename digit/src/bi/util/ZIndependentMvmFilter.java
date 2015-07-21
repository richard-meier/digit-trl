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
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ConcurrentLinkedQueue;

public class ZIndependentMvmFilter implements Runnable{
	private ConcurrentLinkedQueue<StringPair> input;
	private MappingValidation mapVal;
	private double pairedValidityRatioThreshold;
	private int proximityThreshold;
	private GenomeRandomFileAccess genomeAccessor;
	private int RLEN;
	File output;
	
	private long progress;
	
	public static boolean useCustomSimulationRemovalSequence=false;
	public static double[] span;
	
	public ZIndependentMvmFilter(ConcurrentLinkedQueue<StringPair> input, File output, GenomeRandomFileAccess genomeAccessor, double pairedValidityRatioThreshold, int proximityThreshold, int RLEN){
		mapVal = new MappingValidation(pairedValidityRatioThreshold);
		this.pairedValidityRatioThreshold = pairedValidityRatioThreshold;
		this.proximityThreshold = proximityThreshold;
		this.RLEN=RLEN;
		this.genomeAccessor=genomeAccessor;
		this.input=input;
		this.output=output;
		this.progress=0;
	}
	
	@Override
	public void run() {
		try {
			independentMvmFilter(input,output);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void independentMvmFilter(ConcurrentLinkedQueue<StringPair> in, File out) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedWriter allValuesWriter = new BufferedWriter(new FileWriter(out));
		String line1, line2;

		while(in.peek()!=null){
			StringPair pair=in.poll();
			line1=pair.string1;
			line2=pair.string2;
			progress+=1;
			SAMEntry r1 = new SAMEntry(line1);
			SAMEntry r2 = new SAMEntry(line2);
			
			if(!r1.getQname().equals(r2.getQname())){
				System.err.println("ERROR: unpaired reads in thread safe queue!");
				System.err.println("\t"+r1.getQname());
				System.err.println("\t"+r2.getQname());
				System.exit(-1);
			}
			if(doesValidityFilterFail(r1,r2,proximityThreshold,allValuesWriter)){
				if(mapVal.lastValue<0) {
					continue;
				}
				continue;
			}
			if(line1.contains("\tZM:Z:")){
				bw.write(line1);
				bw.newLine();
				bw.write(line2);
				bw.newLine();
			}
			else{
				bw.write(line1+"\tZM:Z:"+mapVal.lastValue);
				bw.newLine();
				bw.write(line2+"\tZM:Z:"+mapVal.lastValue);
				bw.newLine();
			}
		}
		bw.close();
		allValuesWriter.close();
	}
	
	public long getProgress(){
		return progress;
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
			System.err.println("\tSkipping read pair!");
			return true;
		}
		lowerValidity = mapVal.calculateLowerValidityScoreAndPrint(seq1, seq2, reg1, reg2, bw);
		if(useCustomSimulationRemovalSequence){
			if(lowerValidity<span[0] || lowerValidity>span[1]) return false;
			else return true;
		}
		else{
			if(mapVal.checklowerValidityScore(lowerValidity, "custom")) return false;
			else return true;
		}
	}
}
