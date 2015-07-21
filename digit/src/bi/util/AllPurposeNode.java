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

public class AllPurposeNode {
	private AllPurposeNode top;
	private AllPurposeNode left,right,next,previous;
	private double center;
	
	private int position;
	private String chromosomeName;
	private ReadTag tag;
	private AllPurposeNode mate;
	private boolean isRead;
	
	private char strand;
	
	private boolean leadsToReadLeafs;
	
	public boolean traversed=false;
	
	public AllPurposeNode(int position, ReadTag tag, String chromosomeNumber, char strand){
		this.position=position;
		this.center = position;
		this.chromosomeName=chromosomeNumber;
		this.tag=tag;
		left=null;
		right=null;
		isRead=true;
		this.strand=strand;
	}
	
	public AllPurposeNode(int position, ReadTag tag, String chromosomeNumber){
		this.position=position;
		this.center = position;
		this.chromosomeName=chromosomeNumber;
		this.tag=tag;
		left=null;
		right=null;
		isRead=true;
		this.strand='#';
	}
	
	public AllPurposeNode(double center){
		this.center = center;
		left=null;
		right=null;
		isRead=false;
	}
	
	public AllPurposeNode createReadClone(){
		if(!isRead) return null;
		AllPurposeNode out;
		if(this.hasStrand()){
			out = new AllPurposeNode(position, tag, chromosomeName, strand);
		}
		else{
			out = new AllPurposeNode(position, tag, chromosomeName);
		}
		out.setMate(this.mate);
		return out;
	}
	
	public void invokeReadPath(){
		leadsToReadLeafs=true;
		if(top==null) return;
		if(top.leadsToReads()) return;
		else{
			top.invokeReadPath();
		}
	}
	
	private static int testCount;
	public static void connectLeafsAsALinkedList(AllPurposeNode treeRoot, AllPurposeNode listRoot, int allReads){
		testCount=0;
		connect(treeRoot,listRoot, allReads);
	}
	
	private static AllPurposeNode connect(AllPurposeNode subTree, AllPurposeNode lastListElement, int allReads){
		if(subTree==null) {
			return lastListElement;
		}
		if(!subTree.leadsToReads()) {
			return lastListElement;
		}
		testCount++;
		if(allReads>1000){
			if(testCount % (allReads/20) == 0){
				double perc = Math.round(100.0*testCount/allReads);
				System.out.println("\t\t"+perc+"% processed ...");
			}
		}
		if(subTree.isRead()){
			AllPurposeNode previous = connect(subTree.getLeft(),lastListElement, allReads);
			previous.setNext(subTree);
			subTree.setPrevious(previous);
			return connect(subTree.getRight(), subTree, allReads);
		}
		else{
			AllPurposeNode previous = connect(subTree.getLeft(),lastListElement, allReads);
			return connect(subTree.getRight(), previous, allReads);
		}
	}
	
	public boolean leadsToReads(){
		return leadsToReadLeafs;
	}

	public boolean isRead(){
		return isRead;
	}
	
	public boolean hasStrand(){
		return strand == '+' || strand == '-';
	}
	
	public AllPurposeNode getTop(){
		return top;
	}
	
	public void setTop(AllPurposeNode top){
		this.top=top;
	}
	
	public AllPurposeNode getLeft() {
		return left;
	}

	public void setLeft(AllPurposeNode left) {
		this.left = left;
	}

	public AllPurposeNode getRight() {
		return right;
	}

	public void setRight(AllPurposeNode right) {
		this.right = right;
	}
	
	public double getCenter(){
		return center;
	}
	
	public ReadTag getTag(){
		return tag;
	}

	public ReadTag getTagAndAddPositionalInformationToIt(){
		if(tag.getPositionTag()==null){
			tag.setPositionTag(new PositionTag(chromosomeName,position));
		}
		else{}
		return tag;
	}

	public String getChromosomeName() {
		return chromosomeName;
	}

	public AllPurposeNode getMate() {
		return mate;
	}

	public void setMate(AllPurposeNode mate) {
		this.mate = mate;
	}
	
	public AllPurposeNode getPrevious() {
		return previous;
	}

	public void setPrevious(AllPurposeNode previous) {
		this.previous = previous;
	}

	public AllPurposeNode getNext() {
		return next;
	}

	public void setNext(AllPurposeNode next) {
		this.next = next;
	}
	
	public int getPosition(){
		return position;
	}
	
	public char getStrand(){
		return strand;
	}
}
