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

public class ProximityCheckZ {
	private AllPurposeNode orderRoot;
	private AllPurposeNode searchRoot;
	private HashMap<String, AllPurposeNode> readToPosition;
	private int proximityThreshold;
	
	public ProximityCheckZ(int proximityThreshold, int maxPosition){
		orderRoot = new AllPurposeNode(-1,null,null,'/');
		readToPosition = new HashMap<String, AllPurposeNode>();
		this.proximityThreshold = proximityThreshold;
		initialiseSearchTree(maxPosition);
	}
	
	public void printPositionMap(){
		AllPurposeNode curPos = orderRoot;
		while( (curPos=curPos.getNext())!=null ){
			System.out.println(curPos+"\t"+curPos.getTag().getName()+curPos.getTag().getFileNumber()+"\t"+curPos.getPosition()+"\t"+curPos.getChromosomeName());
		}
	}
	
	public ReadTag getTagByName(String name, int number){
		if(number==1) return readToPosition.get(name+"1").getTag();
		else if(number==2) return readToPosition.get(name+"2").getTag();
		else return null;
	}
	
	public AllPurposeNode getAllPurposeNodeByName(String name, int number){
		if(number==1) return readToPosition.get(name+"1");
		else if(number==2) return readToPosition.get(name+"2");
		else return null;
	}
	
	public void setProximityThreshold(int threshold){
		proximityThreshold=threshold;
	}
	
	private void initialiseSearchTree(int maxPos){
		int nFactorial = 1;
		Vector<AllPurposeNode> oldLevel= new Vector<AllPurposeNode>();
		Vector<AllPurposeNode> newLevel= new Vector<AllPurposeNode>();
		AllPurposeNode temp;
		int removeLevelFromBottom = 3;
		while(nFactorial<=maxPos){
			nFactorial=nFactorial*2;
		}
		nFactorial/=2;
		double minLevel = (1.0*maxPos)/nFactorial;
		minLevel*=Math.pow(2, removeLevelFromBottom);
		
		for(double i = minLevel; i<=maxPos; i+=minLevel){
			oldLevel.add(new AllPurposeNode(i-(minLevel/2)));
		}
		while(oldLevel.size()>1){
			for(int i=0;i<oldLevel.size();i+=2){
				temp = new AllPurposeNode((oldLevel.get(i).getCenter()+oldLevel.get(i+1).getCenter())/2);
				temp.setLeft(oldLevel.get(i));
				temp.setRight(oldLevel.get(i+1));
				temp.getLeft().setTop(temp);
				temp.getRight().setTop(temp);
				newLevel.add(temp);
			}
			oldLevel = newLevel;
			newLevel = new Vector<AllPurposeNode>();
		}
		searchRoot = oldLevel.get(0);
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
			AllPurposeNode rankOne = new AllPurposeNode(pos1,r1,chr1,strand1),
						rankTwo = new AllPurposeNode(pos2,r2,chr2,strand2);
			rankOne.setMate(rankTwo);
			rankTwo.setMate(rankOne);
			addB(r1, rankOne);
			addB(r2, rankTwo);
		}
	}
	
	private void addB(ReadTag name, AllPurposeNode entry){
		readToPosition.put(name.getName()+name.getFileNumber(), entry);
		AllPurposeNode start = searchRoot;
		AllPurposeNode nextMember = entry.getPosition()<=start.getLeft().getCenter() ?
			start.getLeft() : start.getRight();
		boolean leftInactive=false;
		boolean rightInactive=false;
		while(nextMember !=null){
			start = nextMember;
			leftInactive=start.getLeft()==null;
			rightInactive=start.getRight()==null;
			if(leftInactive && rightInactive){
				nextMember=null;
				break;
			}
			else if(leftInactive){
				nextMember=start.getRight();
			}
			else if(rightInactive){
				nextMember=start.getLeft();
			}
			else{
				nextMember = entry.getPosition()<=start.getLeft().getCenter() ?
					start.getLeft() : start.getRight();
			}
		}
		if(entry.getPosition()<=start.getCenter()){
			start.setLeft(entry);
		}
		else{
			start.setRight(entry);
		}
		entry.setTop(start);
		entry.invokeReadPath();
	}
	
	public void finishTree(int allReads){
		AllPurposeNode.connectLeafsAsALinkedList(searchRoot, orderRoot, allReads);
	}
	
	public AllPurposeNode getAllPurposeNode(String name){
		return readToPosition.get(name);
	}
	
	public String toString(){
		String out2="D_f: ";
		String out="D_f: ";
		AllPurposeNode newPrev=orderRoot;
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
		AllPurposeNode origin = readToPosition.get(readName+"1");
		AllPurposeNode mate = readToPosition.get(readName+"2");
		mate = origin.getMate();
		emptyList1.add(origin.getTagAndAddPositionalInformationToIt());
		emptyList2.add(mate.getTagAndAddPositionalInformationToIt());
		findProximityCandidates(emptyList1,origin);
		findProximityCandidates(emptyList2,mate);
		trimMateCluster(origin,emptyList1,emptyList2);
	}
	
	private void trimMateCluster(AllPurposeNode originA, Vector<ReadTag> a, Vector<ReadTag> b){
		int oldSize1 = a.size();
		int oldSize2 = b.size();
		do{
			oldSize1 = a.size();
			oldSize2 = b.size();
			refineCluster(originA,a,b);
			refineCluster(originA.getMate(),b,a);
		} while(a.size()!=oldSize1 || b.size()!=oldSize2);
	}
	
	private void refineCluster(AllPurposeNode origin, Vector<ReadTag> originClusterCandidates, Vector<ReadTag> mateClusterCandidates){
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
	
	private void trimClusterRightEnd(AllPurposeNode origin, Vector<ReadTag> originClusterCandidates, int originIndex){
		if(originIndex != originClusterCandidates.size()-1){
			int index = originIndex+1;
			int rightMostBasePosition = origin.getPosition();
			boolean thresholdWasExceeded = false;
			while(index < originClusterCandidates.size()){
				ReadTag currentTag = originClusterCandidates.get(index);
				AllPurposeNode currentRank = readToPosition.get(currentTag.getName()+currentTag.getFileNumber());
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
	
	private void trimClusterLeftEnd(AllPurposeNode origin, Vector<ReadTag> originClusterCandidates, int originIndex){
		if(originIndex != 0){
			int index = originIndex-1;
			int leftMostBasePosition = origin.getPosition();
			boolean thresholdWasExceeded = false;
			while(index >= 0){
				ReadTag currentTag = originClusterCandidates.get(index);
				AllPurposeNode currentRank = readToPosition.get(currentTag.getName()+currentTag.getFileNumber());
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
	
	private void findProximityCandidates(Vector<ReadTag> emptyList, AllPurposeNode readRank0){
		int ankerPoint = readRank0.getPosition();
		AllPurposeNode readRank = readRank0;
		int orPos = readRank.getPosition();
		String orChr = readRank.getChromosomeName();
		AllPurposeNode curPrev = readRank.getPrevious();
		int cpPos = curPrev.getPosition();
		String cpChr = curPrev.getChromosomeName();
		while(cpPos>-1){
			cpPos = curPrev.getPosition();
			cpChr = curPrev.getChromosomeName();
			if(Math.abs(orPos-cpPos)>2*proximityThreshold){
				break;
			}
			if(Math.abs(ankerPoint-cpPos)>10000000){
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
		
		AllPurposeNode curNxt = readRank.getNext();
		if(curNxt == null) return;
		int cnPos;
		String cnChr;
		while(curNxt!=null){
			cnPos = curNxt.getPosition();
			cnChr = curNxt.getChromosomeName();
			if(Math.abs(orPos-cnPos)>2*proximityThreshold){
				break;
			}
			if(Math.abs(ankerPoint-cnPos)>10000000){
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
