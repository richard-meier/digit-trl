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
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import bi.util.AnnotationCheck;
import bi.util.ReadOrientation;
import bi.util.ReadPair;
import bi.util.RegionTag;

public class FicoreLibrary {
	private static final String GROUP_PREFIX="INTERNAL_GROUP_";
	private Vector<String> files;
	private HashMap<String,Vector<ReadPair>> chrPermToReadPair;
	private Vector<String> diseaseFiles,normalFiles;
	private Vector<Vector<String>> groups;
	private Vector<Vector<ReadPair>> countToClusterRegions;
	private String chromosomeFileString;
	Vector<ReadPair> specificOverlapRegions, predispositionRegions;
	AnnotationCheck checkRepeats;
	String outputFile;
	
	private int minimalSampleCount=1;
	
	private static String[] chromosomes;
	private static String[] chromosomes0 = {"hs1","hs10","hs10_GL383545v1_alt","hs10_GL383546v1_alt","hs10_KI270824v1_alt","hs10_KI270825v1_alt","hs11","hs11_GL383547v1_alt","hs11_JH159136v1_alt","hs11_JH159137v1_alt","hs11_KI270721v1_random",
		"hs11_KI270826v1_alt","hs11_KI270827v1_alt","hs11_KI270829v1_alt","hs11_KI270830v1_alt","hs11_KI270831v1_alt","hs11_KI270832v1_alt","hs11_KI270902v1_alt","hs11_KI270903v1_alt","hs11_KI270927v1_alt","hs12","hs12_GL383549v1_alt",
		"hs12_GL383550v2_alt","hs12_GL383551v1_alt","hs12_GL383552v1_alt","hs12_GL383553v2_alt","hs12_GL877875v1_alt","hs12_GL877876v1_alt","hs12_KI270833v1_alt","hs12_KI270834v1_alt","hs12_KI270835v1_alt","hs12_KI270836v1_alt",
		"hs12_KI270837v1_alt","hs12_KI270904v1_alt","hs13","hs13_KI270838v1_alt","hs13_KI270839v1_alt","hs13_KI270840v1_alt","hs13_KI270841v1_alt","hs13_KI270842v1_alt","hs13_KI270843v1_alt","hs14","hs14_GL000009v2_random",
		"hs14_GL000194v1_random","hs14_GL000225v1_random","hs14_KI270722v1_random","hs14_KI270723v1_random","hs14_KI270724v1_random","hs14_KI270725v1_random","hs14_KI270726v1_random","hs14_KI270844v1_alt","hs14_KI270845v1_alt",
		"hs14_KI270846v1_alt","hs14_KI270847v1_alt","hs15","hs15_GL383554v1_alt","hs15_GL383555v2_alt","hs15_KI270727v1_random","hs15_KI270848v1_alt","hs15_KI270849v1_alt","hs15_KI270850v1_alt","hs15_KI270851v1_alt","hs15_KI270852v1_alt",
		"hs15_KI270905v1_alt","hs15_KI270906v1_alt","hs16","hs16_GL383556v1_alt","hs16_GL383557v1_alt","hs16_KI270728v1_random","hs16_KI270853v1_alt","hs16_KI270854v1_alt","hs16_KI270855v1_alt","hs16_KI270856v1_alt","hs17",
		"hs17_GL000205v2_random","hs17_GL000258v2_alt","hs17_GL383563v3_alt","hs17_GL383564v2_alt","hs17_GL383565v1_alt","hs17_GL383566v1_alt","hs17_JH159146v1_alt","hs17_JH159147v1_alt","hs17_JH159148v1_alt","hs17_KI270729v1_random",
		"hs17_KI270730v1_random","hs17_KI270857v1_alt","hs17_KI270858v1_alt","hs17_KI270859v1_alt","hs17_KI270860v1_alt","hs17_KI270861v1_alt","hs17_KI270862v1_alt","hs17_KI270907v1_alt","hs17_KI270908v1_alt","hs17_KI270909v1_alt",
		"hs17_KI270910v1_alt","hs18","hs18_GL383567v1_alt","hs18_GL383568v1_alt","hs18_GL383569v1_alt","hs18_GL383570v1_alt","hs18_GL383571v1_alt","hs18_GL383572v1_alt","hs18_KI270863v1_alt","hs18_KI270864v1_alt","hs18_KI270911v1_alt",
		"hs18_KI270912v1_alt","hs19","hs19_GL000209v2_alt","hs19_GL383573v1_alt","hs19_GL383574v1_alt","hs19_GL383575v2_alt","hs19_GL383576v1_alt","hs19_GL949746v1_alt","hs19_GL949747v2_alt","hs19_GL949748v2_alt","hs19_GL949749v2_alt",
		"hs19_GL949750v2_alt","hs19_GL949751v2_alt","hs19_GL949752v1_alt","hs19_GL949753v2_alt","hs19_KI270865v1_alt","hs19_KI270866v1_alt","hs19_KI270867v1_alt","hs19_KI270868v1_alt","hs19_KI270882v1_alt","hs19_KI270883v1_alt",
		"hs19_KI270884v1_alt","hs19_KI270885v1_alt","hs19_KI270886v1_alt","hs19_KI270887v1_alt","hs19_KI270888v1_alt","hs19_KI270889v1_alt","hs19_KI270890v1_alt","hs19_KI270891v1_alt","hs19_KI270914v1_alt","hs19_KI270915v1_alt",
		"hs19_KI270916v1_alt","hs19_KI270917v1_alt","hs19_KI270918v1_alt","hs19_KI270919v1_alt","hs19_KI270920v1_alt","hs19_KI270921v1_alt","hs19_KI270922v1_alt","hs19_KI270923v1_alt","hs19_KI270929v1_alt","hs19_KI270930v1_alt",
		"hs19_KI270931v1_alt","hs19_KI270932v1_alt","hs19_KI270933v1_alt","hs19_KI270938v1_alt","hs1_GL383518v1_alt","hs1_GL383519v1_alt","hs1_GL383520v2_alt","hs1_KI270706v1_random","hs1_KI270707v1_random","hs1_KI270708v1_random",
		"hs1_KI270709v1_random","hs1_KI270710v1_random","hs1_KI270711v1_random","hs1_KI270712v1_random","hs1_KI270713v1_random","hs1_KI270714v1_random","hs1_KI270759v1_alt","hs1_KI270760v1_alt","hs1_KI270761v1_alt","hs1_KI270762v1_alt",
		"hs1_KI270763v1_alt","hs1_KI270764v1_alt","hs1_KI270765v1_alt","hs1_KI270766v1_alt","hs1_KI270892v1_alt","hs2","hs20","hs20_GL383577v2_alt","hs20_KI270869v1_alt","hs20_KI270870v1_alt","hs20_KI270871v1_alt","hs21","hs21_GL383578v2_alt",
		"hs21_GL383579v2_alt","hs21_GL383580v2_alt","hs21_GL383581v2_alt","hs21_KI270872v1_alt","hs21_KI270873v1_alt","hs21_KI270874v1_alt","hs22","hs22_GL383582v2_alt","hs22_GL383583v2_alt","hs22_KB663609v1_alt","hs22_KI270731v1_random",
		"hs22_KI270732v1_random","hs22_KI270733v1_random","hs22_KI270734v1_random","hs22_KI270735v1_random","hs22_KI270736v1_random","hs22_KI270737v1_random","hs22_KI270738v1_random","hs22_KI270739v1_random","hs22_KI270875v1_alt",
		"hs22_KI270876v1_alt","hs22_KI270877v1_alt","hs22_KI270878v1_alt","hs22_KI270879v1_alt","hs22_KI270928v1_alt","hs2_GL383521v1_alt","hs2_GL383522v1_alt","hs2_GL582966v2_alt","hs2_KI270715v1_random","hs2_KI270716v1_random",
		"hs2_KI270767v1_alt","hs2_KI270768v1_alt","hs2_KI270769v1_alt","hs2_KI270770v1_alt","hs2_KI270771v1_alt","hs2_KI270772v1_alt","hs2_KI270773v1_alt","hs2_KI270774v1_alt","hs2_KI270775v1_alt","hs2_KI270776v1_alt","hs2_KI270893v1_alt",
		"hs2_KI270894v1_alt","hs3","hs3_GL000221v1_random","hs3_GL383526v1_alt","hs3_JH636055v2_alt","hs3_KI270777v1_alt","hs3_KI270778v1_alt","hs3_KI270779v1_alt","hs3_KI270780v1_alt","hs3_KI270781v1_alt","hs3_KI270782v1_alt",
		"hs3_KI270783v1_alt","hs3_KI270784v1_alt","hs3_KI270895v1_alt","hs3_KI270924v1_alt","hs3_KI270934v1_alt","hs3_KI270935v1_alt","hs3_KI270936v1_alt","hs3_KI270937v1_alt","hs4","hs4_GL000008v2_random","hs4_GL000257v2_alt",
		"hs4_GL383527v1_alt","hs4_GL383528v1_alt","hs4_KI270785v1_alt","hs4_KI270786v1_alt","hs4_KI270787v1_alt","hs4_KI270788v1_alt","hs4_KI270789v1_alt","hs4_KI270790v1_alt","hs4_KI270896v1_alt","hs4_KI270925v1_alt","hs5",
		"hs5_GL000208v1_random","hs5_GL339449v2_alt","hs5_GL383530v1_alt","hs5_GL383531v1_alt","hs5_GL383532v1_alt","hs5_GL949742v1_alt","hs5_KI270791v1_alt","hs5_KI270792v1_alt","hs5_KI270793v1_alt","hs5_KI270794v1_alt","hs5_KI270795v1_alt",
		"hs5_KI270796v1_alt","hs5_KI270897v1_alt","hs5_KI270898v1_alt","hs6","hs6_GL000250v2_alt","hs6_GL000251v2_alt","hs6_GL000252v2_alt","hs6_GL000253v2_alt","hs6_GL000254v2_alt","hs6_GL000255v2_alt","hs6_GL000256v2_alt","hs6_GL383533v1_alt",
		"hs6_KB021644v2_alt","hs6_KI270758v1_alt","hs6_KI270797v1_alt","hs6_KI270798v1_alt","hs6_KI270799v1_alt","hs6_KI270800v1_alt","hs6_KI270801v1_alt","hs6_KI270802v1_alt","hs7","hs7_GL383534v2_alt","hs7_KI270803v1_alt","hs7_KI270804v1_alt",
		"hs7_KI270805v1_alt","hs7_KI270806v1_alt","hs7_KI270807v1_alt","hs7_KI270808v1_alt","hs7_KI270809v1_alt","hs7_KI270899v1_alt","hs8","hs8_KI270810v1_alt","hs8_KI270811v1_alt","hs8_KI270812v1_alt","hs8_KI270813v1_alt","hs8_KI270814v1_alt",
		"hs8_KI270815v1_alt","hs8_KI270816v1_alt","hs8_KI270817v1_alt","hs8_KI270818v1_alt","hs8_KI270819v1_alt","hs8_KI270820v1_alt","hs8_KI270821v1_alt","hs8_KI270822v1_alt","hs8_KI270900v1_alt","hs8_KI270901v1_alt","hs8_KI270926v1_alt","hs9",
		"hs9_GL383539v1_alt","hs9_GL383540v1_alt","hs9_GL383541v1_alt","hs9_GL383542v1_alt","hs9_KI270717v1_random","hs9_KI270718v1_random","hs9_KI270719v1_random","hs9_KI270720v1_random","hs9_KI270823v1_alt","hsM","hsUn_GL000195v1",
		"hsUn_GL000213v1","hsUn_GL000214v1","hsUn_GL000216v2","hsUn_GL000218v1","hsUn_GL000219v1","hsUn_GL000220v1","hsUn_GL000224v1","hsUn_GL000226v1","hsUn_KI270302v1","hsUn_KI270303v1","hsUn_KI270304v1","hsUn_KI270305v1","hsUn_KI270310v1",
		"hsUn_KI270311v1","hsUn_KI270312v1","hsUn_KI270315v1","hsUn_KI270316v1","hsUn_KI270317v1","hsUn_KI270320v1","hsUn_KI270322v1","hsUn_KI270329v1","hsUn_KI270330v1","hsUn_KI270333v1","hsUn_KI270334v1","hsUn_KI270335v1","hsUn_KI270336v1",
		"hsUn_KI270337v1","hsUn_KI270338v1","hsUn_KI270340v1","hsUn_KI270362v1","hsUn_KI270363v1","hsUn_KI270364v1","hsUn_KI270366v1","hsUn_KI270371v1","hsUn_KI270372v1","hsUn_KI270373v1","hsUn_KI270374v1","hsUn_KI270375v1","hsUn_KI270376v1",
		"hsUn_KI270378v1","hsUn_KI270379v1","hsUn_KI270381v1","hsUn_KI270382v1","hsUn_KI270383v1","hsUn_KI270384v1","hsUn_KI270385v1","hsUn_KI270386v1","hsUn_KI270387v1","hsUn_KI270388v1","hsUn_KI270389v1","hsUn_KI270390v1","hsUn_KI270391v1",
		"hsUn_KI270392v1","hsUn_KI270393v1","hsUn_KI270394v1","hsUn_KI270395v1","hsUn_KI270396v1","hsUn_KI270411v1","hsUn_KI270412v1","hsUn_KI270414v1","hsUn_KI270417v1","hsUn_KI270418v1","hsUn_KI270419v1","hsUn_KI270420v1","hsUn_KI270422v1",
		"hsUn_KI270423v1","hsUn_KI270424v1","hsUn_KI270425v1","hsUn_KI270429v1","hsUn_KI270435v1","hsUn_KI270438v1","hsUn_KI270442v1","hsUn_KI270448v1","hsUn_KI270465v1","hsUn_KI270466v1","hsUn_KI270467v1","hsUn_KI270468v1","hsUn_KI270507v1",
		"hsUn_KI270508v1","hsUn_KI270509v1","hsUn_KI270510v1","hsUn_KI270511v1","hsUn_KI270512v1","hsUn_KI270515v1","hsUn_KI270516v1","hsUn_KI270517v1","hsUn_KI270518v1","hsUn_KI270519v1","hsUn_KI270521v1","hsUn_KI270522v1","hsUn_KI270528v1",
		"hsUn_KI270529v1","hsUn_KI270530v1","hsUn_KI270538v1","hsUn_KI270539v1","hsUn_KI270544v1","hsUn_KI270548v1","hsUn_KI270579v1","hsUn_KI270580v1","hsUn_KI270581v1","hsUn_KI270582v1","hsUn_KI270583v1","hsUn_KI270584v1","hsUn_KI270587v1",
		"hsUn_KI270588v1","hsUn_KI270589v1","hsUn_KI270590v1","hsUn_KI270591v1","hsUn_KI270593v1","hsUn_KI270741v1","hsUn_KI270742v1","hsUn_KI270743v1","hsUn_KI270744v1","hsUn_KI270745v1","hsUn_KI270746v1","hsUn_KI270747v1","hsUn_KI270748v1",
		"hsUn_KI270749v1","hsUn_KI270750v1","hsUn_KI270751v1","hsUn_KI270752v1","hsUn_KI270753v1","hsUn_KI270754v1","hsUn_KI270755v1","hsUn_KI270756v1","hsUn_KI270757v1","hsX","hsX_KI270880v1_alt","hsX_KI270881v1_alt","hsX_KI270913v1_alt",
		"hsY","hsY_KI270740v1_random"
		};
	
	public static void replaceChromosomeNames(String nameFile) throws IOException{
		Vector<String> names = new Vector<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(nameFile)));
		String line;
		while( (line=br.readLine())!=null ){
			if(line.length()<1) continue;
			names.add(line);
		}
		br.close();
		Collections.sort(
			names, new Comparator<String>(){
				public int compare(String f1, String f2){
					return f1.toString().compareTo(f2.toString());
				}        
			}
		);
		chromosomes=new String[names.size()];
		for(int i=0; i<names.size(); i++){
			chromosomes[i]=names.get(i);
		}
	}

	public FicoreLibrary(String dFileStr, String nFileStr, String chromosomeFile, String[] categoriesStr, String outputFile) throws IOException{
		this.outputFile=outputFile;
		this.chromosomeFileString=chromosomeFile;
//		if(chromosomes!=null){
			replaceChromosomeNames(chromosomeFile);
//		}
		
		files = new Vector<String>();
		normalFiles = new Vector<String>();
		diseaseFiles = new Vector<String>();
		String[] dFiles = dFileStr.split(",");
		String[] nFiles = nFileStr.split(",");
		initialise();
		
		System.err.println("LCc = "+chromosomes.length);
		System.err.println("LC0 = "+chromosomes0.length);
		int foundElements=0;
		for(String chr1 : chromosomes){
			for(String chr2 : chromosomes0){
				if(chr1.equals(chr2)){foundElements++;}
			}
		}
		System.err.println("Number of elements in both lists = "+foundElements);
		
		groups = new Vector<Vector<String>>();
		if(categoriesStr!=null && categoriesStr.length>0) for(String category:categoriesStr){
			String[] files = category.split(",");
			Vector<String> nCat = new Vector<String>();
			groups.add(nCat);
			addFileCategory(files,nCat,false);
		}

		addFileCategory(dFiles,diseaseFiles, true);
		addFileCategory(nFiles,normalFiles, true);
		
		for(int z=0; z<groups.size(); z++){
			Vector<String> group = groups.get(z);
			String groupIdentifier = GROUP_PREFIX+z; 
			for(int i=0; i<group.size(); i++){
				for(int j=0; j<diseaseFiles.size(); j++){
					if(group.get(i).equals(diseaseFiles.get(j))){
						if(!diseaseFiles.contains(groupIdentifier)) diseaseFiles.add(groupIdentifier);
						diseaseFiles.remove(i);
						break;
					}
				}
				for(int j=0; j<normalFiles.size(); j++){
					if(group.get(i).equals(normalFiles.get(j))){
						if(!normalFiles.contains(groupIdentifier)) normalFiles.add(groupIdentifier);
						normalFiles.remove(i);
						break;
					}
				}
			}
		}
		
		System.out.println("ALL TRL COUNT :: "+countR);
		System.out.println("SPECIFIED_FILE_GROUPS: "+groups.size());
		for(int u=0; u<groups.size(); u++){
			System.out.print("\t"+GROUP_PREFIX+u+"\tsize="+groups.get(u).size()+"\t");
			System.out.print(groups.get(u).get(0));
			for(int v=1; v<groups.get(u).size(); v++){
				System.out.print(","+groups.get(u).get(v));
			}
			System.out.println();
		}
		System.out.println("DISEASE_SAMPLES: "+diseaseFiles.size());
		System.out.println("NORMAL_SAMPLES: "+normalFiles.size());
		System.out.println("ALL_SAMPLES: "+files.size());
		System.out.println();
		countToClusterRegions = new Vector<Vector<ReadPair>>();
		for(int cnt=0;cnt<files.size()+2;cnt++){
			countToClusterRegions.add(new Vector<ReadPair>());
		}
		specificOverlapRegions = new Vector<ReadPair>();
		predispositionRegions = new Vector<ReadPair>();
	}
	
	private void addFileCategory(String[] fileNames, Vector<String> category, boolean addToHashMaps) throws IOException{
		for(String s:fileNames){
			if(s.equals("ficore")) continue;
			String fs=s.charAt(0)=='@'?s.substring(1, s.length()):s;
			File current=new File(fs);
			if(!current.exists()){
				System.out.println("FILE DOES NOT EXIST: "+s);
				continue;
			}
			if(addToHashMaps) {
				if(s.charAt(0)=='@'){
					fillHashMapsWithLibrary(chromosomes, new File(fs),category,false);
				}
				else fillHashMapsWithLibrary(chromosomes, new File(fs),category,true);
			}
			else{
				category.add(fs);
			}
		}
	}
	
	private void initialise() throws IOException{
//		BufferedReader br=new BufferedReader(new FileReader("/home/richard/Desktop/new_jobFiles/ficore_dump/hg38_chr.txt")); aqewdsaesffcaqewasfsdq	
		BufferedReader br=new BufferedReader(new FileReader(chromosomeFileString));
		String line;
		Vector<String> chrs=new Vector<String>();
		while((line=br.readLine())!=null){
			if(line.length()<2) continue;
			if(line.charAt(0)=='c' && line.charAt(1)=='h' && line.charAt(2)=='r'){
				line="hs"+line.substring(3);
			}
			else if(line.charAt(0)!='h' && line.charAt(1)!='s'){
				line="hs"+line;
			}
			chrs.add(line);
		}
		Collections.sort(chrs);
		chromosomes=new String[chrs.size()];
		for(int i=0;i<chrs.size();i++) {
			chromosomes[i]=chrs.get(i);
		}
		chrPermToReadPair = new HashMap<String,Vector<ReadPair>>();
		for(int i=0; i<chromosomes.length;i++){
			for(int j=i+1;j<chromosomes.length;j++){
				chrPermToReadPair.put(chromosomes[i]+chromosomes[j], new Vector<ReadPair>());
			}
		}
		br.close();
	}
	
	private void fillHashMapsWithLibrary(String[] chromosomes, File lib, Vector<String> category, boolean addPrefix) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(lib));
		String line;
		while( (line=br.readLine())!=null ){
			if(line.length()==0) continue;
			if(line.charAt(0)=='#') continue;
			line=line.replace("chr", "hs");
			String[] entries = line.split("\\s");
			ReadPair pair = new ReadPair();
			String chrom1,chrom2;
			int start1,start2,stop1,stop2;
			boolean swap=false;
			if(entries[0].compareTo(entries[3])<0){
				chrom1=entries[0];
				start1=Integer.parseInt(entries[1]);
				stop1=Integer.parseInt(entries[2]);
				chrom2=entries[3];
				start2=Integer.parseInt(entries[4]);
				stop2=Integer.parseInt(entries[5]);
			}
			else{
				chrom1=entries[3];
				start1=Integer.parseInt(entries[4]);
				stop1=Integer.parseInt(entries[5]);
				chrom2=entries[0];
				start2=Integer.parseInt(entries[1]);
				stop2=Integer.parseInt(entries[2]);
				swap=true;
			}

			boolean in1=false, in2=false;
			for(String chr:chromosomes){
				if(chr.equals(chrom1)){
					in1=true;
				}
				if(chr.equals(chrom2)){
					in2=true;
				}
			}
			if(!in1 || !in2) continue;

			pair.reg1=new RegionTag(chrom1,start1,stop1);
			pair.reg2=new RegionTag(chrom2,start2,stop2);
			pair.iDs=new Vector<String>();
			pair.fileNameToReadsInRegion1=new HashMap<String,Vector<ReadOrientation>>();
			pair.fileNameToReadsInRegion2=new HashMap<String,Vector<ReadOrientation>>();
			
			String p1 = br.readLine();
			String p2 = br.readLine();
			String prefix="";
			if(addPrefix) {
				prefix=lib.getAbsolutePath()+":";
			}
			
			updateRegionSubFileCount(pair, p1, p2, category, prefix, swap);

			chrPermToReadPair.get(chrom1+chrom2).add(pair);
			countR++;
		}
		br.close();
	}
	
	private static long countR=0;
	
	private void updateRegionSubFileCount(ReadPair region, String p1, String p2, Vector<String> category, String fixedFilePrefix, boolean swap) {
		String[] entries1 = p1.split("\\s");
		String[] entries2 = p2.split("\\s");
		Vector<String> pe1 = new Vector<String>();
		Vector<String> pe2 = new Vector<String>();
		for(int i=2; i<entries1.length;i++){
			String full=fixedFilePrefix+entries1[i];
			for(int u=0; u<groups.size(); u++){
				boolean doBreak=false;
				for(int v=0; v<groups.get(u).size(); v++){
					if(full.contains(groups.get(u).get(v))){
						doBreak=true;
						full=GROUP_PREFIX+u+":"+full;
						break;
					}
				}
				if(doBreak) break;
			}
			pe1.add(full);
		}
		for(int i=2; i<entries2.length;i++){
			String full=fixedFilePrefix+entries2[i];
			for(int u=0; u<groups.size(); u++){
				boolean doBreak=false;
				for(int v=0; v<groups.get(u).size(); v++){
					if(full.contains(groups.get(u).get(v))){
						doBreak=true;
						full=GROUP_PREFIX+u+":"+full;
						break;
					}
				}
				if(doBreak) break;
			}
			pe2.add(full);
		}
		if(swap){
			Vector<String> pe0=pe1;
			pe1=pe2;
			pe2=pe0;
		}
		Collections.sort(pe1); Collections.sort(pe2);

		for(int i=0;i<pe1.size();i++){
			String[] s1 = pe1.get(i).split(":");
			ReadOrientation read = new ReadOrientation(); read.reg=region.reg1; read.pos=Integer.parseInt(s1[s1.length-1]); read.orientation=s1[s1.length-3]; read.complete=pe1.get(i);
			String iD=s1[0];
			if(!files.contains(iD)) {
				files.add(iD);
				category.add(iD);
			}
			boolean notInside=!region.fileNameToReadsInRegion1.keySet().contains(iD);
			if(notInside) {
				if(!region.iDs.contains(s1[0])) region.iDs.add(s1[0]);
				region.fileNameToReadsInRegion1.put(iD, new Vector<ReadOrientation>());
			}
			region.fileNameToReadsInRegion1.get(iD).add(read);
		}
		for(int i=0;i<pe2.size();i++){
			String[] s2 = pe2.get(i).split(":");
			ReadOrientation read = new ReadOrientation(); read.reg=region.reg2; read.pos=Integer.parseInt(s2[s2.length-1]); read.orientation=s2[s2.length-3]; read.complete=pe2.get(i);
			String iD=s2[0];
			boolean notInside=!region.fileNameToReadsInRegion2.keySet().contains(iD);
			if(notInside) {
				if(!region.iDs.contains(s2[0])) region.iDs.add(s2[0]);
				region.fileNameToReadsInRegion2.put(iD, new Vector<ReadOrientation>());
			}
			region.fileNameToReadsInRegion2.get(iD).add(read);
		}
	}
	
	public void run() throws IOException{
		findCommonRegions();
	}
	
	private void findCommonRegions() throws IOException{
		int count=0;
		int count2=0;
		int count3=0;
		for(int i=0; i<chromosomes.length;i++){
			for(int j=i+1;j<chromosomes.length;j++){
				Vector<ReadPair> regions = new Vector<ReadPair>();
				count++;
				if(chrPermToReadPair.get(chromosomes[i]+chromosomes[j]).size()==0) continue;
				for(ReadPair rp:chrPermToReadPair.get(chromosomes[i]+chromosomes[j])){
					boolean covered=false;
					for(ReadPair region:regions){
						if(rp.overlaps(region, 10000)){
							covered=true;
							merge(region,rp);
						}
					}
					if(!covered){
						regions.add(rp);
					}
				}
				
				tryToFinalMerge(regions);
				for(ReadPair reg:regions){
					if(reg.iDs.size()>=minimalSampleCount){
						storeCommonRegion(reg);
						count2++;
					}
					count3++;
				}
			}
		}
System.err.println(count+"\t"+chrPermToReadPair.keySet().size()+"\t"+count2+"\t"+count3);
		System.out.println("==========< ALL UNIQUE >==========");
		reportAllUnique();
		System.out.println();
		System.out.println("============< COMMON >============");
		reportAllCommon();
		System.out.println();
		System.out.println("========< PREDISPOSITION >========");
		reportPredispositions(false);
		System.out.println();
		
		System.out.println("===========< SPECIFIC >===========");
		reportAllCancerSpecific(false);
		System.out.println();
		System.out.println("=========< SPECIFIC-SUB >=========");
		reportAllCancerSpecific(true);
//		System.out.println();
		System.out.println("=========< SUPER TARGET >=========");
		findCommonSuperClusterTargets();
		System.out.println();
		if(outputFile!=null){
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputFile)));
			generateNewLibrary(bw);
			bw.close();
		}
		System.out.println("==========< BREAK DOWN >==========");
		System.out.println("Common = " + allCommon);
		System.out.println("Specific = "+allSpecific);
		
		reportClusterPropertiesOnly();
	}
	
	private static int allCommon=0, allSpecific=0;
	
	private void merge(ReadPair target, ReadPair input){
		int start1 = Math.min(target.reg1.getStart(), input.reg1.getStart());
		int start2 = Math.min(target.reg2.getStart(), input.reg2.getStart());
		int stop1 = Math.max(target.reg1.getStop(), input.reg1.getStop());
		int stop2 = Math.max(target.reg2.getStop(), input.reg2.getStop());
		target.reg1.setStart(start1);
		target.reg1.setStop(stop1);
		target.reg2.setStart(start2);
		target.reg2.setStop(stop2);
		for(String fStr:input.iDs){
			if(!target.iDs.contains(fStr)) target.iDs.add(fStr);
			if(!target.fileNameToReadsInRegion1.keySet().contains(fStr)){
				target.fileNameToReadsInRegion1.put(fStr, new Vector<ReadOrientation>());
			}
			for(ReadOrientation read:input.fileNameToReadsInRegion1.get(fStr)){
				target.fileNameToReadsInRegion1.get(fStr).add(read);
			}
			if(!target.fileNameToReadsInRegion2.keySet().contains(fStr)){
				target.fileNameToReadsInRegion2.put(fStr, new Vector<ReadOrientation>());
			}
			for(ReadOrientation read:input.fileNameToReadsInRegion2.get(fStr)){
				target.fileNameToReadsInRegion2.get(fStr).add(read);
			}
		}
	}
	
	private void tryToFinalMerge(Vector<ReadPair> regions){
		for(int i=0; i<regions.size(); i++){
			for(int j=i+1; j<regions.size(); j++){
				if(regions.get(i).overlaps(regions.get(j), 10000)){
					merge(regions.get(i),regions.get(j));
				}
			}
		}
	}
	
	private void storeCommonRegion(ReadPair reg) {
		Vector<ReadPair> currentRegion=countToClusterRegions.get(reg.iDs.size());
		currentRegion.add(reg);
	}
	
	private void reportAllUnique(){
		Vector<ReadPair> currentRegion=countToClusterRegions.get(1);
		
		int cancerCount=0;
		for(int j=0; j<currentRegion.size(); j++){
			reportCommonRegion(currentRegion.get(j),"");
			ReadPair reg = currentRegion.get(j);
			int ac = howManyFilesAreShared(reg.iDs,diseaseFiles);
			if(ac>0){
				cancerCount++;
			}
		}
		System.out.println("Total="+currentRegion.size());
		System.out.println("Target="+cancerCount);
	}
	
	private void reportAllCommon(){
		for(int i=countToClusterRegions.size()-1; i>0; i--){
			Vector<ReadPair> currentRegion=countToClusterRegions.get(i);
			int cnt=0;
			for(int j=0; j<currentRegion.size(); j++){
				ReadPair reg = currentRegion.get(j);
				int ac = howManyFilesAreShared(reg.iDs,diseaseFiles);
				int bc = howManyFilesAreShared(reg.iDs,normalFiles);

				if(bc>0 && ac+bc>1){
					reportCommonRegion(currentRegion.get(j),"");
					allCommon++;
					if(i>0) specificOverlapRegions.add(currentRegion.get(j));
					cnt++;
				}

			}
			if(cnt>0) System.out.println();
		}
	}
	
	private void reportPredispositions(boolean reportSub) throws IOException{
		Vector<SortTupel> tups = new Vector<SortTupel>();
		for(int i=countToClusterRegions.size()-1; i>0; i--){
			Vector<ReadPair> currentRegion=countToClusterRegions.get(i);
			for(int j=0; j<currentRegion.size(); j++){
				ReadPair reg = currentRegion.get(j);

				int cancerCount = howManyFilesAreShared(reg.iDs,diseaseFiles);
				int normalCount = howManyFilesAreShared(reg.iDs,normalFiles);
				if(normalCount==0 || cancerCount==0) continue;

				double normalizedCancerCount=1.0*cancerCount/(1.0*diseaseFiles.size());
				double normalizedNormalCount=1.0*normalCount/(1.0*normalFiles.size());
				double cancerOccuranceRate = Math.log10(normalizedCancerCount/normalizedNormalCount);
				if(cancerOccuranceRate<0.204119 /*0.19477*/ /*0.1760912*/) continue;
				SortTupel tup = new SortTupel();
				tup.pair=currentRegion.get(j);
				tup.rate=cancerOccuranceRate;
				tups.add(tup);
				if(i>0) predispositionRegions.add(currentRegion.get(j));
			}
		}
		Collections.sort(tups);
		for(int i=0; i<tups.size(); i++){
			SortTupel tup = tups.get(i);
			double rounded = Math.round(tup.rate*10000.0)/10000.0;
			reportCommonRegion(tup.pair, rounded+"\t");
			reportCommonRegionWithProperties(tup.pair);
		}
	}
	
	private void reportAllCancerSpecific(boolean reportSub) throws IOException{
		for(int i=countToClusterRegions.size()-1; i>0; i--){
			Vector<ReadPair> currentRegion=countToClusterRegions.get(i);
			int cnt=0;
			for(int j=0; j<currentRegion.size(); j++){
				ReadPair reg = currentRegion.get(j);
				
				int ac = howManyFilesAreShared(reg.iDs,diseaseFiles);
				int bc = howManyFilesAreShared(reg.iDs,normalFiles);

				if( ac>0 && bc==0){
					if(reportSub) {
						reportCommonRegionWithProperties(currentRegion.get(j));
					}
					else {
						reportCommonRegion(currentRegion.get(j),"");
					}
					allSpecific++;
					if(i>0) specificOverlapRegions.add(currentRegion.get(j));
					cnt++;
				}

			}
			if(cnt>0) System.out.println();
		}
	}
	
	private void reportClusterPropertiesOnly() throws IOException{
		for(int i=countToClusterRegions.size()-1; i>0; i--){
			Vector<ReadPair> currentRegion=countToClusterRegions.get(i);
			for(int j=0; j<currentRegion.size(); j++){
				reportSubClusterProperties(currentRegion.get(j));
			}
		}
	}
	
	private void reportCommonRegion(ReadPair reg, String prefix) {
		System.out.println(prefix+reg.reg1.getChromosome()+":"+reg.reg1.getStart()+"-"+reg.reg1.getStop()+"\t"+
			reg.reg2.getChromosome()+":"+reg.reg2.getStart()+"-"+reg.reg2.getStop()+"\t"+reg.iDs.size()+"\t"+reg.iDs);
	}
	
	private void reportCommonRegionWithProperties(ReadPair reg) throws IOException{
		HashMap<String,Integer> orientationCounter = new HashMap<String,Integer>();
		orientationCounter.put("++", 0);
		orientationCounter.put("+-", 0);
		orientationCounter.put("-+", 0);
		orientationCounter.put("--", 0);
		String infoOut="";
		for(String s:reg.fileNameToReadsInRegion1.keySet()){
			Vector<ReadOrientation> rds1 = reg.fileNameToReadsInRegion1.get(s);
			Vector<ReadOrientation> rds2 = reg.fileNameToReadsInRegion2.get(s);
			int start1=Integer.MAX_VALUE,start2=Integer.MAX_VALUE;
			int stop1=Integer.MIN_VALUE,stop2=Integer.MIN_VALUE;
			String orientations="";
			for(int i=0; i<rds1.size();i++){
				ReadOrientation rd1=rds1.get(i);
				ReadOrientation rd2=rds2.get(i);
				start1=rd1.pos<start1?rd1.pos:start1;
				start2=rd2.pos<start2?rd2.pos:start2;
				stop1=rd1.pos>stop1?rd1.pos:stop1;
				stop2=rd2.pos>stop2?rd2.pos:stop2;
				String orientation=rd1.orientation+rd2.orientation;

				orientationCounter.put(orientation, orientationCounter.get(orientation)+1);
				orientations+=orientation+" ";
			}
			infoOut+="\t"+reg.reg1.getChromosome()+":"+start1+"-"+stop1+"\t"+reg.reg2.getChromosome()+":"+start2+"-"+stop2+"\t"+s+"\t"+orientations+"\n";
		}
		System.out.println(reg.reg1.getChromosome()+":"+reg.reg1.getStart()+"-"+reg.reg1.getStop()+"\t"+
			reg.reg2.getChromosome()+":"+reg.reg2.getStart()+"-"+reg.reg2.getStop()+"\t"+reg.iDs.size()+"\t"+
			"(++,--,+-,-+)=("+orientationCounter.get("++")+","+orientationCounter.get("--")+","+orientationCounter.get("+-")+","+orientationCounter.get("-+")+")"+"\n"+infoOut);
	}
	
	private void reportSubClusterProperties(ReadPair reg) throws IOException{
		HashMap<String,Integer> orientationCounter = new HashMap<String,Integer>();
		orientationCounter.put("++", 0);
		orientationCounter.put("+-", 0);
		orientationCounter.put("-+", 0);
		orientationCounter.put("--", 0);
		for(String s:reg.fileNameToReadsInRegion1.keySet()){
			Vector<ReadOrientation> rds1 = reg.fileNameToReadsInRegion1.get(s);
			Vector<ReadOrientation> rds2 = reg.fileNameToReadsInRegion2.get(s);
			int start1=Integer.MAX_VALUE,start2=Integer.MAX_VALUE;
			int stop1=Integer.MIN_VALUE,stop2=Integer.MIN_VALUE;
			for(int i=0; i<rds1.size();i++){
				ReadOrientation rd1=rds1.get(i);
				ReadOrientation rd2=rds2.get(i);
				start1=rd1.pos<start1?rd1.pos:start1;
				start2=rd2.pos<start2?rd2.pos:start2;
				stop1=rd1.pos>stop1?rd1.pos:stop1;
				stop2=rd2.pos>stop2?rd2.pos:stop2;
				String orientation=rd1.orientation+rd2.orientation;

				orientationCounter.put(orientation, orientationCounter.get(orientation)+1);
			}
		}
	}
	
	private int howManyFilesAreShared(Vector<String> list1, Vector<String> list2){
		int count=0;
		for(String f1:list1){
			for(String f2: list2){
				if(f1.equals(f2)) count++;
			}
		}
		return count;
	}
	
	private void findCommonSuperClusterTargets(){
		Vector<ReadPair> allSuperClusters = new Vector<ReadPair>();
		for(int i=0; i<countToClusterRegions.size(); i++){
			for(ReadPair reg:countToClusterRegions.get(i)) {
				int cancerCount = howManyFilesAreShared(reg.iDs,diseaseFiles);
				int normalCount = howManyFilesAreShared(reg.iDs,normalFiles);
				if(cancerCount==0) continue;
				if(normalCount>0){
					double normalizedCancerCount=1.0*cancerCount/(1.0*diseaseFiles.size());
					double normalizedNormalCount=1.0*normalCount/(1.0*normalFiles.size());
					double cancerOccuranceRate = Math.log10(normalizedCancerCount/normalizedNormalCount);
					if(cancerOccuranceRate<0.204119 /*0.19477*/ /*0.1760912*/) continue;
				}
				allSuperClusters.add(reg);
			}
		}
		Vector<SuperTarget> allTargets = new Vector<SuperTarget>();
		Vector<ReadPair> doneClusters = new Vector<ReadPair>();
		for(int i=0; i<allSuperClusters.size(); i++){
			ReadPair cluster = allSuperClusters.get(i);
			if(doneClusters.contains(cluster)){
				continue;
			}
			SuperTarget st1=new SuperTarget(cluster.reg1,cluster);
			SuperTarget st2=new SuperTarget(cluster.reg2,cluster);
			for(int j=0; j<allSuperClusters.size(); j++){
				if(i==j) continue;
				if(doneClusters.contains(cluster)){
					continue;
				}
				ReadPair currentOrigin=allSuperClusters.get(j);
				if(currentOrigin==cluster) continue;
				RegionTag t1=currentOrigin.reg1;
				RegionTag t2=currentOrigin.reg2;
				if(
						st1.tryToAddAndUpdate(t1, currentOrigin) ||
						st1.tryToAddAndUpdate(t2, currentOrigin) ||
						st2.tryToAddAndUpdate(t1, currentOrigin) ||
						st2.tryToAddAndUpdate(t2, currentOrigin)
				){
					doneClusters.add(currentOrigin);
				}
			}
			if(st1.targets.size()>1){
				allTargets.add(st1);
			}
			if(st2.targets.size()>1){
				allTargets.add(st2);
			}
			doneClusters.add(cluster);
		}
		for(int i=35; i>0; i--){
			boolean done=false;
			for(SuperTarget target:allTargets){
				int count=0;
				for(ReadPair rp:target.superClusterOrigins){
					count+=rp.iDs.size();
				}
				if(count!=i) continue;
				done=true;
				System.out.println(target.chr+":"+target.start+"-"+target.stop+"\t"+count);
				for(ReadPair rp:target.superClusterOrigins){
					System.out.print("\t");
					reportCommonRegion(rp,"");
				}
			}
			if(done) System.out.println();
		}
	}
	
	private void generateNewLibrary(BufferedWriter bw) throws IOException{
		for(int i=countToClusterRegions.size()-1; i>0; i--){
			Vector<ReadPair> currentRegion=countToClusterRegions.get(i);
			int cnt=0;
			for(int j=0; j<currentRegion.size(); j++){
				ReadPair reg = currentRegion.get(j);
				int ac = howManyFilesAreShared(reg.iDs,diseaseFiles);
				int bc = howManyFilesAreShared(reg.iDs,normalFiles);
				if( ac>=0 && bc>0 && ac+bc>=2){
					printLibraryOutPut(currentRegion.get(j), bw);
					cnt++;
				}
			}
			if(cnt>0) bw.newLine();
		}
	}
	
	private void printLibraryOutPut(ReadPair reg, BufferedWriter bw) throws IOException{
		Vector<ReadOrientation> allReads1=new Vector<ReadOrientation>();
		Vector<ReadOrientation> allReads2=new Vector<ReadOrientation>();
		for(String s:reg.fileNameToReadsInRegion1.keySet()){
			Vector<ReadOrientation> rds1 = reg.fileNameToReadsInRegion1.get(s);
			Vector<ReadOrientation> rds2 = reg.fileNameToReadsInRegion2.get(s);
			for(int i=0; i<rds1.size();i++){
				ReadOrientation rd1=rds1.get(i);
				ReadOrientation rd2=rds2.get(i);
				allReads1.add(rd1);
				allReads2.add(rd2);
			}
		}
		String outMain = 	reg.reg1.getChromosome()+" "+reg.reg1.getStart()+" "+reg.reg1.getStop()+" "+
				reg.reg2.getChromosome()+" "+reg.reg2.getStart()+" "+reg.reg2.getStop();
		String outP1 = "#\t"+allReads1.size(), 
				outP2 = "#\t"+allReads2.size();
		for(ReadOrientation rd:allReads1){
			outP1+="\t"+rd.complete;
		}
		for(ReadOrientation rd:allReads2){
			outP2+="\t"+rd.complete;
		}
		bw.write(outMain); bw.newLine();
		bw.write(outP1); bw.newLine();
		bw.write(outP2); bw.newLine();
	}
	
	
	private class SortTupel implements Comparable<SortTupel>{
		public ReadPair pair;
		public double rate;
		@Override
		public int compareTo(SortTupel other) {
			if(this.rate>other.rate) return -1;
			else if(this.rate<other.rate) return 1;
			else return 0;
		}
	}
	
	private class SuperTarget{
		public Vector<RegionTag> targets;
		public Vector<ReadPair> superClusterOrigins;
		public String chr;
		public int start, stop;
		
		public SuperTarget(RegionTag target, ReadPair origin){
			targets=new Vector<RegionTag>();
			targets.add(target);
			superClusterOrigins=new Vector<ReadPair>();
			superClusterOrigins.add(origin);
			start = target.getStart();
			stop = target.getStop();
			chr = target.getChromosome();
		}
		
		public boolean tryToAddAndUpdate(RegionTag target, ReadPair origin){
			if(target.getChromosome() == null) return false;
			if(target.getChromosome().equals(chr)){
				if(overlaps(target,11000)){
					start=Math.min(target.getStart(), start);
					stop=Math.max(target.getStop(), stop);
					targets.add(target);
					superClusterOrigins.add(origin);
					return true;
				}
			}
			return false;
		}
		
		public boolean overlaps(RegionTag other, int outerRange){
			boolean overlap = 
				(start>= other.getStart()-outerRange && start<= other.getStop()+outerRange) ||
				(stop>= other.getStart()-outerRange && stop<= other.getStop()+outerRange) ||
				(start<= other.getStart()-outerRange && stop>= other.getStop()+outerRange) ||
				(start>= other.getStart()-outerRange && stop<= other.getStop()+outerRange);
			return overlap;
		}
	}
	
	
}
