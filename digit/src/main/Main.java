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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import bi.util.GenomeRandomFileAccess;

public class Main {

	public static String 
	inputFile = null, 
	inputFileR2 = null,
	inputFile2 = null,
	inputFile3 = null,
	inputFile4 = null,
	outputPath = null,
	outputPath2 = null,
	chromLengthFile = null;
	
	public static String[] inputArray1=null;
	public static double numberTimesSigma;
	public static int proximityCount;
	public static int proximityThreshold, threadNumber;
	private static int readLength;
	private static int anker;
	private static double mvmConcordantThreshold,mapQual;
	private static String type;
	public static boolean printReadNames;
	public static double[] enteredMvmSpan;
	public static boolean fqIsOn=true;
	
	public static final String u_index = "index\n"+
			"       -g genome     \t specify genome fasta file\n"+
			"       -o output     \t specify output file for the created index\n\n";
	
	public static final String u_ficore = "ficore\n"+
			"       -D  target    \t specify list of files in target (disease) group\n"+
			"       -N  other     \t specify list of files in other (normal) group\n"+
			"       -C  groups    \t specify sub-groups\n"+
			"       -r  chr-file  \t chromosome name file for reference\n"+
			"       -lc output   \t specify output file for library creation\n\n"
			;
	
	public static final String u_analyse = "analyse\n"+
			"       -i  infile   \t specify SAM-input file\n"+
			"       -o  outpath  \t specify directory for output files\n"+
			"       -s  sigma    \t specify the multiple of sigma [2.33]\n"+
			"       -lc lc-file  \t bed file with low complexity region annotation of the genome\n"+
			"       -q  quality  \t MAPQ threshold used to reject reads [0]\n"+
			"       -r  chr-file \t chromosome length file for reference\n"+
			"       -a  aln      \t filter by XT:A:U tags (legacy setting) [false]\n\n";

	public static final String u_revcomp = "revcomp\n"+
			"       -i infile \t specify fastq input file\n"+
			"       -o outfile \t specify fastq output file\n\n";

	public static final String u_proxval = "proxval\n"+
			"       -i infile     \t specify analysed input file\n"+
			"       -o outfile    \t specify output file\n"+
			"       -c count      \t specify the proximity count used for validation [3] \n"+
			"       -t threshold  \t specify the proximity threshold used for validation [100]\n"+
			"       -g genome     \t specify genome fasta file\n"+
			"       -x index      \t specify genome index file\n"+
			"       -r chr-file   \t chromosome length file for reference\n"+
			"       -M MVM        \t specify the MVM concordant threshold [0.005]\n"+
			"       -T rd-length  \t specify the maximum read length\n\n";
	
	public static final String u_mvmfo = "mvmfo\n"+
			"       -i infile     \t specify analysed input file\n"+
			"       -o outfile    \t specify output file\n"+
			"       -t threshold  \t specify the proximity threshold used for validation [100]\n"+
			"       -g genome     \t specify genome fasta file\n"+
			"       -x index      \t specify genome index file\n"+
			"       -M MVM        \t specify the MVM concordant threshold [0.005]\n"+
			"       -p thread-num \t specify the number of threads to be used\n"+
			"       -r chr-file   \t chromosome length file for reference\n"+
			"       -T rd-length  \t specify the maximum read length\n\n"
			;

	public static final String u_chrospl = "chrospl\n"+
			"       -i  input     \t specify headerless SAM-input file\n"+
			"       -o1 output1   \t specify primary output file\n"+
			"       -o2 output2   \t specify secondary output file\n"+
			"       -s  chrset    \t specify chromosome set that input is split by\n\n"
			;
	
	public static final String usage=
			"\n"+
			u_analyse+
			u_chrospl+
			u_index+
			u_mvmfo+
			u_proxval+
			u_revcomp
			;

	public static void main(String[] args) throws IOException {
		numberTimesSigma=2.33;
		if(args.length==0){
			System.out.println(usage);
			System.exit(0);
		}
		
		if(args[0].equals("index")){
			applyIndexInputArguments(args);
			GenomeRandomFileAccess.buildIndex(inputFile, outputPath);
		}
		else if(args[0].equals("revcomp")){
			applyRevCompInputArguments(args);
			ReverseComplementer revComp = new ReverseComplementer(new File(inputFile), new File(outputPath), fqIsOn);
			revComp.run();
		}
		else if(args[0].equals("analyse")){
			mapQual=0;
			chromLengthFile=null;
			applyAnalyseInputArguments(args);
			Analyser myNewAnalyser=new Analyser(inputFile, outputPath, numberTimesSigma, type, new File(inputFile2), mapQual, chromLengthFile);
			myNewAnalyser.run();
		}
		else if(args[0].equals("proxval")){
			anker=-1;
			readLength=-1;
			proximityCount = 3;
			proximityThreshold = 100;
			printReadNames=false;
			chromLengthFile=null;
			applyProxvalArguments(args);
			ProximityValidation proxVal=new ProximityValidation(new File(inputFile), new File(outputPath), proximityCount, proximityThreshold, inputFile2, inputFile3, true, mvmConcordantThreshold, chromLengthFile);
			if(anker>0){
				proxVal.overwriteAnkerLimit(anker);
				System.out.println("USING alternative anchor span: "+anker);
			}
			proxVal.performClusterAnalysis(readLength);
		}
		else if(args[0].equals("mvmfo")){
			readLength=-1;
			proximityCount = 3;
			proximityThreshold = 100;
			printReadNames=false;
			mvmConcordantThreshold=-1;
			threadNumber=1;
			chromLengthFile=null;
			applyPureFilterArguments(args);
			ProximityValidation proxVal=new ProximityValidation(new File(inputFile), new File(outputPath), proximityCount, proximityThreshold, inputFile2, inputFile3, false, mvmConcordantThreshold, chromLengthFile);
			try {
				if(enteredMvmSpan !=null){
					proxVal.performMvmFilter(new File(inputFile),new File(outputPath),readLength,enteredMvmSpan,threadNumber);
				}
				else proxVal.performMvmFilter(new File(inputFile),new File(outputPath),readLength,null,threadNumber);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		else if(args[0].equals("ficore")){
			readLength=500;
			applyFicoreArguments(args);
			FicoreLibrary fcr = new FicoreLibrary(inputFile,inputFile2,inputFile3,inputArray1,outputPath);
			fcr.run();
		}
		else if(args[0].equals("chrospl")){
			applyChromosomeSplittingArguments(args);
			ChromosomeSplitting csp = new ChromosomeSplitting(inputArray1,inputFile,outputPath,outputPath2);
			csp.run();
		}
		else if(args[0].equals("dev.test.sw")){
			long tseed = Long.parseLong(args[1]);
			String storeOut = null;
			if(args.length>2){
				if(args[2].length()>1){
					storeOut = args[2];
				} else{
					System.err.println("ERROR::Invalid file destination!");
					System.exit(-1);
				}
			}
			TestSwAlignment.run(storeOut, tseed);
		}
		else {
			System.err.println("Error:: Invalid command.");
			System.out.println(usage);
			System.exit(-1);
		}
	}

	private static void applyFicoreArguments(String[] args) {
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-lc")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-D")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-N")){
				i++;
				inputFile2=args[i];
			}
			else if(args[i].equals("-C")){
				i++;
				inputArray1=args[i].split(":");
			}
			else if(args[i].equals("-r")){
				i++;
				inputFile3=args[i];
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_ficore);
				System.exit(-1);
			}
		}
		if(inputFile == null || inputFile2 == null) {
			System.err.println("Error:INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_ficore);
			System.exit(-1);
		}
	}

	private static void applyIndexInputArguments(String[] args) {
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-g")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-o")){
				i++;
				outputPath=args[i];
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_index);
				System.exit(-1);
			}
		}
		if(inputFile == null || outputPath == null) {
			System.err.println("Error:INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_index);
			System.exit(-1);
		}
	}

	private static void applyAnalyseInputArguments(String[] args) {
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-i")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-o")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-s")){
				i++;
				numberTimesSigma=Double.parseDouble(args[i]);
			}
			else if(args[i].equals("-q")){
				i++;
				mapQual=Double.parseDouble(args[i]);
			}
			else if(args[i].equals("-a")){
				i++;
				type=args[i];
			}
			else if(args[i].equals("-lc")){
				i++;
				inputFile2=args[i];
			}
			else if(args[i].equals("-r")){
				i++;
				chromLengthFile=args[i];
			}
			else if(args[i].equals("-help")){
				System.err.println(u_analyse);
				System.exit(-1);
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_analyse);
				System.out.println(args[i]);
				System.exit(-1);
			}
		}
		if(inputFile == null || outputPath == null) {
			System.err.println("Error:INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_analyse);
			System.exit(-1);
		}
	}

	private static void applyRevCompInputArguments(String[] args) {
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-i")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-o")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-fa")){
				fqIsOn=false;
			}
			else if(args[i].equals("-help")){
				System.err.println(u_revcomp);
				System.exit(-1);
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_revcomp);
				System.exit(-1);
			}
		}
		if(inputFile == null || outputPath == null) {
			System.err.println("Error: INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_revcomp);
			System.exit(-1);
		}
	}
	
	private static void applyProxvalArguments(String[] args) throws IOException{
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-i")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-g")){
				i++;
				inputFile2=args[i];
			}
			else if(args[i].equals("-x")){
				i++;
				inputFile3=args[i];
			}
			else if(args[i].equals("-o")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-c")){
				i++;
				proximityCount=Integer.parseInt(args[i]);
			}
			else if(args[i].equals("-t")){
				i++;
				if(args[i].contains("file=")){
					String fileString = args[i].replace("file=", "");
					proximityThreshold=parseOutProximityThreshold(new File(fileString));
				}
				else {
					proximityThreshold=Integer.parseInt(args[i]);
				}
			}
			else if(args[i].equals("-r")){
				i++;
				chromLengthFile=args[i];
			}
			else if(args[i].equals("-T")){
				i++;
				readLength=Integer.parseInt(args[i]);
			}
			else if(args[i].equals("-M")){
				i++;
				mvmConcordantThreshold=Double.parseDouble(args[i]);
			}
			else if(args[i].equals("-an")){
				i++;
				anker=Integer.parseInt(args[i]);
			}
			else if(args[i].equals("-help")){
				System.err.println(u_proxval);
				System.exit(-1);
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_proxval);
				System.exit(-1);
			}
		}
		if(inputFile == null || outputPath == null) {
			System.err.println("Error: INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_proxval);
			System.exit(-1);
		}
	}
	
	private static void applyPureFilterArguments(String[] args) throws IOException{
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-i")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-g")){
				i++;
				inputFile2=args[i];
			}
			else if(args[i].equals("-x")){
				i++;
				inputFile3=args[i];
			}
			else if(args[i].equals("-o")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-t")){
				i++;
				if(args[i].contains("file=")){
					String fileString = args[i].replace("file=", "");
					proximityThreshold=parseOutProximityThreshold(new File(fileString));
				}
				else {
					proximityThreshold=Integer.parseInt(args[i]);
				}
			}
			else if(args[i].equals("-T")){
				i++;
				readLength=Integer.parseInt(args[i]);
			}
			else if(args[i].equals("-M")){
				i++;
				mvmConcordantThreshold=Double.parseDouble(args[i]);
			}
			else if(args[i].equals("-sp")){
				i++;
				String[] tmp=args[i].split("_");
				enteredMvmSpan=new double[]{Double.parseDouble(tmp[0]),Double.parseDouble(tmp[1])};
			}
			else if(args[i].equals("-p")){
				i++;
				threadNumber=Integer.parseInt(args[i]);
			}
			else if(args[i].equals("-r")){
				i++;
				chromLengthFile=args[i];
			}
			else if(args[i].equals("-help")){
				System.err.println(u_mvmfo);
				System.exit(-1);
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_mvmfo);
				System.exit(-1);
			}
		}
		
		if(inputFile == null || outputPath == null) {
			System.err.println("Error: INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_mvmfo);
			System.exit(-1);
		}
	}
	
	private static int parseOutProximityThreshold(File file) throws IOException{
		int out=-1;
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while( (line=br.readLine())!=null ){
			if(line.length()<2) {
				continue;
			}
			if(line.contains("recommended cluster threshold: ")){
				out=Integer.parseInt(line.replaceAll("recommended cluster threshold: ", ""));
				break;
			}
		}
		br.close();
		return out;
	}
	
	private static void applyChromosomeSplittingArguments(String[] args) {
		for(int i=1;i<args.length;i++){
			if(args[i].equals("-i")){
				i++;
				inputFile=args[i];
			}
			else if(args[i].equals("-o1")){
				i++;
				outputPath=args[i];
			}
			else if(args[i].equals("-o2")){
				i++;
				outputPath2=args[i];
			}
			else if(args[i].equals("-s")){
				i++;
				inputArray1=args[i].split(",");
			}
			else{
				System.err.println("Error: invalid input arguments, see usage");
				System.out.println(u_chrospl);
				System.exit(-1);
			}
		}
		if(inputFile == null || outputPath == null) {
			System.err.println("Error:INVALID INPUT: Please define ALL non-optional parameters!");
			System.out.println(u_chrospl);
			System.exit(-1);
		}
	}
	
}
