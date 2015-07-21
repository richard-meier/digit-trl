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

import java.util.HashMap;
import java.util.List;
import java.util.Vector;

public class ProximityCheckB {
	private PositionRank currentListBuild;
	private PositionRank orderRoot;
	private PositionFinder searchRoot;
	private HashMap<String, PositionRank> readToPosition;
	private int proximityThreshold;
	public int ankerThreshold=10000000;
	public long readNumberTotal=0;
	
	public ProximityCheckB(int proximityThreshold, int maxPosition){
		orderRoot = new PositionRankSD(-1,null,null,'/');
		readToPosition = new HashMap<String, PositionRank>();
		this.proximityThreshold = proximityThreshold;
		initialiseSearchTree(maxPosition);
		currentListBuild=orderRoot;
	}
	
	public void printPositionMap(){
		PositionRank curPos = orderRoot;
		while( (curPos=curPos.getNext())!=null ){
			System.out.println(curPos+"\t"+curPos.getTag().getName()+curPos.getTag().getFileNumber()+"\t"+curPos.getPosition()+"\t"+curPos.getChromosomeName());
		}
	}
	
	public ReadTag getTagByName(String name, int number){
		if(number==1) return readToPosition.get(name+"1").getTag();
		else if(number==2) return readToPosition.get(name+"2").getTag();
		else return null;
	}
	
	public PositionRank getPositionRankByName(String name, int number){
		if(number==1) return readToPosition.get(name+"1");
		else if(number==2) return readToPosition.get(name+"2");
		else return null;
	}
	
	public void setProximityThreshold(int threshold){
		proximityThreshold=threshold;
	}

	private void initialiseSearchTree(int maxPos){
		int nFactorial = 1;
		Vector<PositionFinder> oldLevel= new Vector<PositionFinder>();
		Vector<PositionFinder> newLevel= new Vector<PositionFinder>();
		PositionFinder temp;
		int removeLevelFromBottom = 3;
		while(nFactorial<=maxPos){
			nFactorial=nFactorial*2;
		}
		nFactorial/=2;
		double minLevel = (1.0*maxPos)/nFactorial;
		minLevel*=Math.pow(2, removeLevelFromBottom);
		
		for(double i = minLevel; i<=maxPos; i+=minLevel){
			oldLevel.add(new PositionFinder(i-(minLevel/2)));
		}
		while(oldLevel.size()>1){
			for(int i=0;i<oldLevel.size();i+=2){
				temp = new PositionFinder((oldLevel.get(i).getCenter()+oldLevel.get(i+1).getCenter())/2);
				temp.setLeft(oldLevel.get(i));
				temp.setRight(oldLevel.get(i+1));
				newLevel.add(temp);
			}
			oldLevel = newLevel;
			newLevel = new Vector<PositionFinder>();
		}
		searchRoot = oldLevel.get(0);
	}
	
	public void addLinear(ReadTag r1, int pos1, String chr1, char strand1){
		if(readToPosition.containsKey(r1)){
			System.err.println("WARNING::Submitted read name is occuring multiple times!\nEntry '"+
				r1+"' will be ignored!");
		}
		else {
			PositionRank rankOne = new PositionRankSD(pos1,r1,chr1,strand1);
			currentListBuild.setNext(rankOne);
			rankOne.setPrevious(currentListBuild);
			currentListBuild=rankOne;
			readToPosition.put(r1.getName()+r1.getFileNumber(), rankOne);
			readNumberTotal++;
		}
	}
	
	public void linkMates(){
		System.out.println("\tLinking Read Mates ...");
		PositionRank currentListElement = orderRoot;
		int pcnt=0;
		while(currentListElement.getNext()!=null){
			currentListElement=currentListElement.getNext();
			if(currentListElement.getMate()==null){
				ReadTag tag = currentListElement.getTag();
				String name2=(tag.getFileNumber()==(byte)1) ? tag.getName()+2 : tag.getName()+1;
				PositionRank mate=readToPosition.get(name2);
				mate.setMate(currentListElement);
				currentListElement.setMate(mate);
			}
			pcnt++;
			if(pcnt%1000000==0){
				System.out.println("\t\t"+pcnt+" reads processed!");
			}
		}
		System.out.println("\tLinking Done!");
	}
	
	public void add(ReadTag r1, ReadTag r2, int pos1, int pos2, String chr1, String chr2, char strand1, char strand2){
		if(readToPosition.containsKey(r1)||readToPosition.containsKey(r2)){
			System.err.println("WARNING::Submitted read names are occuring multiple times!\nEntries '"+
				r1+"' and '"+r2+"' will be ignored!");
		}
		else 
		if(r1.equals(r2)){
			System.err.println("WARNING::Submitted read names have the same value!\nEntries '"+
					r1+"'(1) and '"+r2+"'(2) will be ignored!");
		}
		else{
			PositionRank rankOne = new PositionRankSD(pos1,r1,chr1,strand1),
						rankTwo = new PositionRankSD(pos2,r2,chr2,strand2);
			rankOne.setMate(rankTwo);
			rankTwo.setMate(rankOne);
			addB(r1, rankOne);
			addB(r2, rankTwo);
		}
	}
	
	private void addB(ReadTag name, PositionRank entry){
		readToPosition.put(name.getName()+name.getFileNumber(), entry);
		PositionFinder start = searchRoot;
		PositionRank neighbour = null;
		while((neighbour=start.getEnd())==null && start.getLeft()!=null){
			if(entry.getPosition()<=start.getLeft().getCenter()){
				start=start.getLeft();
			}
			else start=start.getRight();
		}
		if(neighbour != null){
			while(neighbour.getPosition()>entry.getPosition()) neighbour=neighbour.getPrevious();

			while(neighbour.getNext()!=null){
				neighbour=neighbour.getNext();
				if(neighbour.getPosition()>entry.getPosition()) break;
			}
			if(!(neighbour.getPosition()>entry.getPosition())){
				neighbour.setNext(entry);
				entry.setPrevious(neighbour);
			}
			else{
				entry.setNext(neighbour);
				entry.setPrevious(neighbour.getPrevious());
				neighbour.getPrevious().setNext(entry);
				neighbour.setPrevious(entry);
			}
		}
		else{
			start.setEnd(entry);
			add(name,entry);
		}
	}
	
	private void add(ReadTag name, PositionRank entry){
		readToPosition.put(name.getName()+name.getFileNumber(), entry);
		if(orderRoot.getNext() == null){
			orderRoot.setNext(entry);
			entry.setPrevious(orderRoot);
		}
		else{
			PositionRank newNext=orderRoot;
			while(newNext.getNext()!=null){
				newNext=newNext.getNext();
				if(newNext.getPosition()>entry.getPosition()) break;
			}
			if(!(newNext.getPosition()>entry.getPosition())){
				newNext.setNext(entry);
				entry.setPrevious(newNext);
			}
			else{
				entry.setNext(newNext);
				entry.setPrevious(newNext.getPrevious());
				newNext.getPrevious().setNext(entry);
				newNext.setPrevious(entry);
			}
		}
	}
	
	public PositionRank getPositionRank(String name){
		return readToPosition.get(name);
	}
	
	public String toString(){
		String out2="D_f: ";
		String out="D_f: ";
		PositionRank newPrev=orderRoot;
		while(newPrev.getNext()!=null){
			newPrev=newPrev.getNext();
			out+=newPrev.getPosition()+" ";
			out2+=newPrev.getTag().getName()+"."+readToPosition.get(newPrev.getTag().getName()+newPrev.getTag().getFileNumber()).getPosition()+" ";
		}
		out2+="\nD_r: ";
		out+="\nD_r: ";
		while(newPrev.getPrevious()!=null){
			out+=newPrev.getPosition()+" ";
			out2+=newPrev.getTag().getName()+"."+readToPosition.get(newPrev.getTag().getName()+newPrev.getTag().getFileNumber()).getPosition()+" ";;
			newPrev=newPrev.getPrevious();
		}
		out+="\n";
		out2+="\n";
		return out+out2;
	}
		
	public void findClusterNeighbours(Vector<ReadTag> emptyList1, Vector<ReadTag> emptyList2, String readName){
		PositionRank origin = readToPosition.get(readName+"1");
		PositionRank mate = readToPosition.get(readName+"2");
		mate = origin.getMate();
		emptyList1.add(origin.getTagAndAddPositionalInformationToIt());
		emptyList2.add(mate.getTagAndAddPositionalInformationToIt());
		findProximityCandidates(emptyList1,origin);
		findProximityCandidates(emptyList2,mate);
		trimMateCluster(origin,emptyList1,emptyList2);
	}
	
	private void trimMateCluster(PositionRank originA, Vector<ReadTag> a, Vector<ReadTag> b){
		int oldSize1 = a.size();
		int oldSize2 = b.size();
		do{
			oldSize1 = a.size();
			oldSize2 = b.size();
			refineCluster(originA,a,b);
			refineCluster(originA.getMate(),b,a);
		} while(a.size()!=oldSize1 || b.size()!=oldSize2);
	}
	
	private void refineCluster(PositionRank origin, Vector<ReadTag> originClusterCandidates, Vector<ReadTag> mateClusterCandidates){
		retainNamesInFrom(originClusterCandidates, mateClusterCandidates);
		int originIndex=0;
		for(originIndex=0; originIndex<originClusterCandidates.size(); originIndex++){
			if(originClusterCandidates.get(originIndex).equals(origin.getTag())){
				break;
			}
		}
		trimClusterRightEnd(origin, originClusterCandidates, originIndex);
		trimClusterLeftEnd(origin, originClusterCandidates, originIndex);
	}
	
	private void trimClusterRightEnd(PositionRank origin, Vector<ReadTag> originClusterCandidates, int originIndex){
		if(originIndex != originClusterCandidates.size()-1){
			int index = originIndex+1;
			int rightMostBasePosition = origin.getPosition();
			boolean thresholdWasExceeded = false;
			while(index < originClusterCandidates.size()){
				ReadTag currentTag = originClusterCandidates.get(index);
				PositionRank currentRank = readToPosition.get(currentTag.getName()+currentTag.getFileNumber());
				int currentBasePosition = currentRank.getPosition();
				boolean rightMostCandidateExceedsProximityThreshold = Math.abs(rightMostBasePosition-currentBasePosition) > 2*proximityThreshold;
				if( rightMostCandidateExceedsProximityThreshold ){
					thresholdWasExceeded = true;
					break;
				}
				else{
					rightMostBasePosition = currentBasePosition;
				}
				index++;
			}
			if(thresholdWasExceeded) removeElementsInRegion(originClusterCandidates, index, originClusterCandidates.size());
		}
	}
	
	private void removeElementsInRegion(List<ReadTag> list, int start, int stop){
		list.subList(start, stop).clear();
	}
	
	private void trimClusterLeftEnd(PositionRank origin, Vector<ReadTag> originClusterCandidates, int originIndex){
		if(originIndex != 0){
			int index = originIndex-1;
			int leftMostBasePosition = origin.getPosition();
			boolean thresholdWasExceeded = false;
			while(index >= 0){
				ReadTag currentTag = originClusterCandidates.get(index);
				PositionRank currentRank = readToPosition.get(currentTag.getName()+currentTag.getFileNumber());
				int currentBasePosition = currentRank.getPosition();
				boolean leftMostCandidateExceedsProximityThreshold = Math.abs(leftMostBasePosition-currentBasePosition) > 2*proximityThreshold;
				if( leftMostCandidateExceedsProximityThreshold ){
					thresholdWasExceeded = true;
					break;
				}
				else{
					leftMostBasePosition = currentBasePosition;
				}
				index--;
			}
			if(thresholdWasExceeded) removeElementsInRegion(originClusterCandidates, 0, index+1);
		}
	}
	
	private void retainNamesInFrom(Vector<ReadTag> a, Vector<ReadTag> b){
		for(int i=0; i<a.size();i++){
			boolean found=false;
			ReadTag r1 = a.get(i);
			for(ReadTag r2:b){
				if(r1.getName().equals(r2.getName()) && r1.getFileNumber() != r2.getFileNumber()){ //TODO: ZOMFG HEREEEEE
					found=true;
					break;
				}
			}
			if(!found) a.set(i,null);
		}
		boolean check=true;
		while(check){
			check=false;
			for(int i=0; i<a.size();i++){
				if(a.get(i)==null) {
					a.remove(i);
					check=true;
					break;
				}
			}
		}
	}
	
	private void findProximityCandidates(Vector<ReadTag> emptyList, PositionRank readRank0){
		int ankerPoint = readRank0.getPosition();
		PositionRank readRank = readRank0;
		int orPos = readRank.getPosition();
		String orChr = readRank.getChromosomeName();
		PositionRank curPrev = readRank.getPrevious();
		int cpPos = curPrev.getPosition();
		String cpChr = curPrev.getChromosomeName();
		while(cpPos>-1){
			cpPos = curPrev.getPosition();
			cpChr = curPrev.getChromosomeName();
			if(Math.abs(orPos-cpPos)>2*proximityThreshold){
				break;
			}
			if(Math.abs(ankerPoint-cpPos)>ankerThreshold){
				break;
			}
			if(!orChr.equals(cpChr)){
				curPrev = curPrev.getPrevious();
				continue;
			}
			emptyList.add(0, curPrev.getTagAndAddPositionalInformationToIt());
			String targetName = curPrev.getTag().getName()+curPrev.getTag().getFileNumber();
			readRank = readToPosition.get(targetName);
			orPos = readRank.getPosition();
			orChr = readRank.getChromosomeName();
			curPrev = curPrev.getPrevious();
		}
		readRank = readRank0;
		orPos = readRank.getPosition();
		orChr = readRank.getChromosomeName();
		
		PositionRank curNxt = readRank.getNext();
		if(curNxt == null) return;
		int cnPos;
		String cnChr;
		while(curNxt!=null){
			cnPos = curNxt.getPosition();
			cnChr = curNxt.getChromosomeName();
			if(Math.abs(orPos-cnPos)>2*proximityThreshold){
				break;
			}
			if(Math.abs(ankerPoint-cnPos)>ankerThreshold){
				break;
			}
			if(!orChr.equals(cnChr)){
				curNxt = curNxt.getNext();
				continue;
			}
			emptyList.add(curNxt.getTagAndAddPositionalInformationToIt());
			String targetName = curNxt.getTag().getName()+curNxt.getTag().getFileNumber();
			readRank = readToPosition.get(targetName);
			orPos = readRank.getPosition();
			orChr = readRank.getChromosomeName();
			curNxt = curNxt.getNext();
		}
	}
}
