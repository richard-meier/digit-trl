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

package bi.sim;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;

import bi.util.DnaMod;
import bi.util.GenomeRandomFileAccess;
import bi.util.ReadPair;
import bi.util.RegionTag;

public class SimulateMatePairSample {
	static Random rand=new Random();
	public static int totalReadNumber=80000000;
	public static int numberOfBreaks=1000;
	public static Random otherRandom=new Random();
	public static double genomeLength = 3.088269832E9;
	public static Random allRand=new Random();
	private static final int CHR_LENGTH_MIN=4600000;
	private static final int avgSepDist=2300;
	private static final double sd=700;
	static String name="sim_batch08";
	static String outPath=".../simulation/sim_batch08/";
	public static double mutationRate=0.001;
	public static double mut2Rate=0.001;
	private static long mutationCount=0;
	public static int RD_MAX=100, RDMIN=30;
	public static void run() {
		ReadPair rA=new ReadPair();
		rA.reg1=new RegionTag("chr7", 49680251, 49686294);
		rA.reg2=new RegionTag("chr11", 95436223, 95442271);
		specRegs=new ReadPair[]{
				rA
		};
		System.out.println(name);
		try{
			GenomeRandomFileAccess genome = new GenomeRandomFileAccess(".../GRCh38/hg38.fa", ".../GRCh38/hg38.iVIGS.idx", null);
			genome.loadIndex();

			System.out.println("PASS! :)");
			Vector<String> chromosomes = new Vector<String>();
			Vector<Double> chrChance = new Vector<Double>();
			Vector<Double> chrCumulativeChance = new Vector<Double>();
			double sum=0;
			for(String s:genome.chromosomeLengthMap.keySet()){
				if(genome.chromosomeLengthMap.get(s)<CHR_LENGTH_MIN) continue;
				if(s.contains("_alt")) continue;
				chromosomes.add(s);
				chrChance.add(genome.chromosomeLengthMap.get(s)+0.0);
				sum+=genome.chromosomeLengthMap.get(s)+0.0;
			}
			for(int i=0; i<chrChance.size(); i++){
				chrChance.set(i, chrChance.get(i)/sum);
				chrCumulativeChance.add(0.0);
			}
			chrCumulativeChance.set(0,chrChance.get(0));
			for(int i=1; i<chrCumulativeChance.size(); i++){
				chrCumulativeChance.set(i, chrCumulativeChance.get(i-1)+chrChance.get(i));
			}
			chrCumulativeChance.set(chrCumulativeChance.size()-1,1.0);
			Vector<SVT> translocations = new Vector<SVT>();
			getTranslocations(genome,translocations,chromosomes);
			
			BufferedWriter bw1=new BufferedWriter(new FileWriter(new File(outPath+name+"_1.fa")));
			BufferedWriter bw2=new BufferedWriter(new FileWriter(new File(outPath+name+"_2.fa")));
			int diff=500000;
			System.out.println("START GENERATING");
			long id=1;
			long trls=0;
			for(int i=0;i<totalReadNumber;i++){
				String chr = getRandomChr(chrCumulativeChance,chromosomes);
				if(genome.chromosomeLengthMap.get(chr)<CHR_LENGTH_MIN){
					i--;
					continue;
				}
				int pos = rand.nextInt(genome.chromosomeLengthMap.get(chr)-diff)+diff/2;
				int sep = getSepWithThreshold();
				String s1;
				String s2;
				boolean found=false;
				boolean chanceToHitTheChromosome=rand.nextBoolean();
				for(SVT trl:translocations){
					if(
							(chr.equals(trl.chromosome1) && pos+100<=trl.breakPoint1 && pos+sep>=trl.breakPoint1) ||
							(chr.equals(trl.chromosome1) && pos+100<=trl.breakPoint1 && pos+sep>=trl.breakPoint1) ||
							(chr.equals(trl.chromosome2) && pos+100<=trl.breakPoint2 && pos+sep>=trl.breakPoint2) ||
							(chr.equals(trl.chromosome2) && pos+100<=trl.breakPoint2 && pos+sep>=trl.breakPoint2)
					){
						if(chanceToHitTheChromosome){
							trl.readPairs++;
							found=true;
							trls++;
						}
					}
				}

				if(found) continue;
				if(i%10000==0) System.out.println("ReadPairs:"+i+"\tMutations:"+mutationCount+"\tTRLS:"+trls);//"\tFakeTranslocations:"+fakeRegions+"\tCarriedOnFTs:"+repls+"\tlastScore:"+tmpScore);
				int l1 = RD_MAX-otherRandom.nextInt(RDMIN);
				int l2 = RD_MAX-otherRandom.nextInt(RDMIN);
					s1=genome.getGenomicSequence(chr, pos, pos+l1).toUpperCase();
					s2=DnaMod.reverseComplement(genome.getGenomicSequence(chr, pos+sep, pos+sep+l2).toUpperCase());
					StringBuffer bs1=new StringBuffer();
					StringBuffer bs2=new StringBuffer();
					boolean preCheck1=false;
					boolean preCheck2=false;
					
				for(int k=0; k<s1.length(); k++){
					if(preCheck1 && otherRandom.nextDouble()<mut2Rate){
						bs1.append(getRandomNucBase(s1.charAt(k)));
					}
					else if(otherRandom.nextDouble()<mutationRate){
						bs1.append(getRandomNucBase(s1.charAt(k)));
						mutationCount++;
					}
					else{
						bs1.append(s1.charAt(k));
					}
				}
				for(int k=0; k<s2.length(); k++){
					if(preCheck2 && otherRandom.nextDouble()<mut2Rate){
						bs2.append(getRandomNucBase(s2.charAt(k)));
					}
					else if(otherRandom.nextDouble()<mutationRate){
						bs2.append(getRandomNucBase(s2.charAt(k)));
						mutationCount++;
					}
					else{
						bs2.append(s2.charAt(k));
					}
				}
				String header=">SIMULATED_READ_PAIR#"+id+":"+chr+":"+(pos)+":"+chr+":"+(pos+sep);
				bw1.write(header); bw1.newLine();
				bw1.write(bs1.toString()); bw1.newLine();
				bw2.write(header); bw2.newLine();
				bw2.write(bs2.toString()); bw2.newLine();
				id++;
			}
			BufferedWriter bwt = new BufferedWriter(new FileWriter(new File(outPath+"breaks_min_N.txt")));
			BufferedWriter bwt2 = new BufferedWriter(new FileWriter(new File(outPath+"breaks_circos.txt")));
			
			for(SVT trl:translocations){
				boolean print=true;
				Vector<String[]> rds = new Vector<String[]>();
				Vector<String> sl1 = new Vector<String>();
				Vector<String> sl2 = new Vector<String>();
				trl.readPairs=(int)(Math.round(trl.readPairs*0.5));
				for(int i=0; i<trl.readPairs; i++){
					int sep=getSepWithThreshold();
					int l1 = RD_MAX-otherRandom.nextInt(RDMIN);
					int l2 = RD_MAX-otherRandom.nextInt(RDMIN);
					int shift=105+allRand.nextInt(sep-211);
					String s1=genome.getGenomicSequence(trl.chromosome1, trl.breakPoint1-shift, trl.breakPoint1-shift+l1).toUpperCase();
					String s2=genome.getGenomicSequence(trl.chromosome2, trl.breakPoint2-shift+sep, trl.breakPoint2-shift+sep+l2).toUpperCase();
					int oInd=rand.nextInt(trl.orientation1.length);
					if(trl.orientation1[oInd].equals("-")){
						s1=DnaMod.reverseComplement(s1);
					}
					if(trl.orientation2[oInd].equals("-")){
						s2=DnaMod.reverseComplement(s2);
					}
					if(countNs(s1)>5 || countNs(s2)>5){
						print=false;
						break;
					}
					
					StringBuffer bs1=new StringBuffer();
					StringBuffer bs2=new StringBuffer();
					for(int k=0; k<s1.length(); k++){
						if(otherRandom.nextDouble()<mutationRate){
							bs1.append(getRandomNucBase(s1.charAt(k)));
							mutationCount++;
						}
						else{
							bs1.append(s1.charAt(k));
						}
					}
					for(int k=0; k<s2.length(); k++){
						if(otherRandom.nextDouble()<mutationRate){
							bs2.append(getRandomNucBase(s2.charAt(k)));
							mutationCount++;
						}
						else{
							bs2.append(s2.charAt(k));
						}
					}
					String header=">SIMULATED_READ_PAIR#"+id+":"+trl.chromosome1+":"+(trl.breakPoint1-shift)+":"+trl.chromosome2+":"+(trl.breakPoint2-shift+sep);
					sl1.add(header); sl1.add(bs1.toString());
					sl2.add(header); sl2.add(bs2.toString());
					rds.add(new String[]{header.substring(1),(trl.breakPoint1-shift)+"",(trl.breakPoint2-shift+sep)+"",trl.orientation1[oInd]+"",trl.orientation2[oInd]+""});
					id++;
				}
				if(print) {
					for(int i=0; i<sl1.size();i++){
						bw1.write(sl1.get(i)); bw1.newLine();
						bw2.write(sl2.get(i)); bw2.newLine();
					}
					bwt.write(trl.toString()); bwt.newLine();
					if(trl.readPairs>0){
						bwt2.write(trl.printHs(rds)); bwt2.newLine();
					}
				}
			}
			bwt.close();
			bwt2.close();
			bw1.close();
			bw2.close();
		}catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
		
	public static double tmpScore=-1;
	public static boolean sequencesAreSufficientlySimilar(String s1, String s2){
		double threshold=90*0.5;
		int score=matchScore(s1,s2);
		if(score>threshold){
			tmpScore=score;
			return true;
		}
		else return false;
	}
	
	public static int matchScore(String s1, String s2){
		int maxCount=0;
		String s1b="";
		for(int u=0; u<110; u++){
			s1b+="#";
		}
		s1b=s1b+s1;
		for(int i=0; i<s1b.length(); i++){
			int tempCount=0;
			for(int j=0; j<s2.length(); j++){
				
				if(i+j>=s1b.length()){
					break;
				}
				if(s1b.charAt(i+j)==s2.charAt(j)){
					tempCount++;
				}
			}
			if(tempCount>maxCount) maxCount=tempCount;
		}
		return maxCount;
	}
	
	public static ReadPair[] specRegs;

	public static int countNs(String s){
		int count=0;
		for(int i=0; i<s.length(); i++){
			if(s.charAt(i)=='N') count++;
		}
		return count;
	}
	
	public static char getRandomNucBase(char in){
		char out;
		if(otherRandom.nextBoolean()){
			if(otherRandom.nextBoolean()){
				out= 'A';
			} else {
				out= 'C';
			}
		} else {
			if(otherRandom.nextBoolean()){
				out= 'G';
			} else {
				out= 'T';
			}
		}
		if(in==out) return getRandomNucBase(in);
		else return out;
	}
	
	private static int getSepWithThreshold(){
		double sep = (Math.round(rand.nextGaussian()*sd+avgSepDist));
		if(sep<310) return getSepWithThreshold();
		else return ((int)(sep));
	}
	
	private static void getTranslocations(GenomeRandomFileAccess genome, Vector<SVT> translocations, Vector<String> chromosomes) throws IOException{
		Random rand=new Random();
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPath+"breaks.txt")));
		int diff=500000;
		for(int i=0; i<numberOfBreaks; i++){
			int index = rand.nextInt(chromosomes.size());
			String chr = chromosomes.get(index);
			if(genome.chromosomeLengthMap.get(chr)<CHR_LENGTH_MIN){
				i--;
				continue;
			}
			int pos = rand.nextInt(genome.chromosomeLengthMap.get(chr)-diff)+diff/2;
			SVT trl = new SVT();
			trl.chromosome1=chr;
			trl.breakPoint1=pos;
			int index2;
			String chr2;
			while( true ){
				index2 = rand.nextInt(chromosomes.size());
				chr2 = chromosomes.get(index2);
				if(chr2.equals(chr)){
					continue;
				}
				if(genome.chromosomeLengthMap.get(chr2)<CHR_LENGTH_MIN){
					continue;
				}
				break;
			}
			int pos2 = rand.nextInt(genome.chromosomeLengthMap.get(chr2)-diff)+diff/2;
			trl.chromosome2=chr2;
			trl.breakPoint2=pos2;
			initialiseOrientations(trl);
			translocations.add(trl);
			bw.write(trl.toString()); bw.newLine();
		}
		bw.close();
	}
	
	public static void initialiseOrientations(SVT trl){
		if(allRand.nextBoolean()){
			trl.orientation1=new String[2];
			trl.orientation2=new String[2];
			if(allRand.nextBoolean()){
				trl.orientation1[0]="+";
				trl.orientation1[1]="-";
				trl.orientation2[0]="+";
				trl.orientation2[1]="-";
			} else{
				trl.orientation1[0]="+";
				trl.orientation1[1]="-";
				trl.orientation2[0]="-";
				trl.orientation2[1]="+";
			}
		} else{
			trl.orientation1=new String[1];
			trl.orientation2=new String[1];
			if(allRand.nextBoolean()){
				if(allRand.nextBoolean()){
					trl.orientation1[0]="+";
					trl.orientation2[0]="+";
				} else{
					trl.orientation1[0]="-";
					trl.orientation2[0]="-";
				}
			} else{
				if(allRand.nextBoolean()){
					trl.orientation1[0]="+";
					trl.orientation2[0]="-";
				} else{
					trl.orientation1[0]="-";
					trl.orientation2[0]="+";
				}
			}
		}
	}

	
	private static String getRandomChr(Vector<Double> chrCumulativeChance, Vector<String> chromosomes){
		double area=allRand.nextDouble();
		int i=0;
		for(i=0; i<chrCumulativeChance.size(); i++){
			if(chrCumulativeChance.get(i)>area){
				break;
			}
		}
		return chromosomes.get(i);
	}
	
	public static class SVT{
		public int readPairs=0;
		public int breakPoint1,breakPoint2;
		public String chromosome1,chromosome2;
		public String[] orientation1, orientation2;
		
		public String toString(){
			String out="";
			if(chromosome1.compareTo(chromosome2)<0){
				out+=chromosome1+":"+breakPoint1+" "+chromosome2+":"+breakPoint2+" "+readPairs+" ";
				out+=orientation1[0]+orientation2[0];
				for(int i=1; i<orientation1.length;i++){
					out+=":"+orientation1[i]+orientation2[i];
				}
			}
			else{
				out+=chromosome2+":"+breakPoint2+" "+chromosome1+":"+breakPoint1+" "+readPairs+" ";
				out+=orientation1[0]+orientation2[0];
				for(int i=1; i<orientation2.length;i++){
					out+=":"+orientation2[i]+orientation1[i];
				}
			}
			return out;
		}
		
		public String printHs(Vector<String[]> headers){
			String out="";
			String l1="";
			String l2="";
			if(chromosome1.compareTo(chromosome2)<0){
				out+=chromosome1+" "+breakPoint1+" "+breakPoint1+" "+chromosome2+" "+breakPoint2+" "+breakPoint2;
				l1+="\n#\t"+headers.size()+"\t"+headers.get(0)[0]+"1:"+headers.get(0)[3]+"::"+headers.get(0)[1];
				l2+="\n#\t"+headers.size()+"\t"+headers.get(0)[0]+"2:"+headers.get(0)[4]+"::"+headers.get(0)[2];
				for(int i=1; i<headers.size();i++){
					l1+="\t"+headers.get(i)[0]+"1:"+headers.get(0)[3]+"::"+headers.get(i)[1];
					l2+="\t"+headers.get(i)[0]+"2:"+headers.get(0)[4]+"::"+headers.get(i)[2];
				}
			}
			else{
				out+=chromosome2+" "+breakPoint2+" "+breakPoint2+" "+chromosome1+" "+breakPoint1+" "+breakPoint1;
				l2+="\n#\t"+headers.size()+"\t"+headers.get(0)[0]+"1:"+headers.get(0)[3]+"::"+headers.get(0)[1];
				l1+="\n#\t"+headers.size()+"\t"+headers.get(0)[0]+"2:"+headers.get(0)[4]+"::"+headers.get(0)[2];
				for(int i=1; i<headers.size();i++){
					l2+="\t"+headers.get(i)[0]+"1:"+headers.get(0)[3]+"::"+headers.get(i)[1];
					l1+="\t"+headers.get(i)[0]+"2:"+headers.get(0)[4]+"::"+headers.get(i)[2];
				}
			}
			out=out.replace("chr", "hs");
			return out+l1+l2;
		}
	}
}
