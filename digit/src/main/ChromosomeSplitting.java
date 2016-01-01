package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ChromosomeSplitting {
	File input, out1, out2;
	String[] primarySet;
	private static final long[] magnitudes = { 
		100l, 1000l, 10000l, 100000l, 1000000l, 10000000l, 100000000l, 1000000000l, 10000000000l, 100000000000l, 1000000000000l, 
		10000000000000l, 100000000000000l, 1000000000000000l, 10000000000000000l, 100000000000000000l,1000000000000000000l
	};
	
	public ChromosomeSplitting(String[] set, String in, String out1, String out2){
		this.input=new File(in);
		this.out1=new File(out1);
		this.out2=new File(out2);
		this.primarySet=set;
	}

	public void run() throws IOException{
		System.out.println("input = "+input.getAbsolutePath());
		System.out.println("Starting to split file...");
		long pairCount=0;
		int magnitudeSurpassIndex=0;
		long startingTime = System.currentTimeMillis();
		try {
			BufferedReader br = new BufferedReader(new FileReader(input));
			BufferedWriter s1 = new BufferedWriter(new FileWriter(out1));
			BufferedWriter s2 = new BufferedWriter(new FileWriter(out2));
			String line1,line2;
			while( (line1=br.readLine()) != null){
				if(line1.length()<2){
					continue;
				}
				line2 = br.readLine();
				String[] entries1 = line1.split("\\s");
				String[] entries2 = line2.split("\\s");
				String chr1=entries1[2];
				String chr2=entries2[2];
				boolean pairIsInPrimary=false;
				for(String primary : primarySet){
					if(chr1.equals(primary) || chr2.equals(primary)){
						pairIsInPrimary=true;
						break;
					}
				}
				if(pairIsInPrimary){
					s1.write(line1); s1.newLine();
					s1.write(line2); s1.newLine();
				}
				else{
					s2.write(line1); s2.newLine();
					s2.write(line2); s2.newLine();
				}
				pairCount++;
				if(pairCount > magnitudes[magnitudeSurpassIndex]){
					System.out.println("\t"+magnitudes[magnitudeSurpassIndex]+" read pairs processed [time = "+(System.currentTimeMillis()-startingTime)+"ms] ...");
					magnitudeSurpassIndex++;
				}
			}
			s1.close();
			s2.close();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("\t"+pairCount+" read pairs processed [time = "+(System.currentTimeMillis()-startingTime)+"ms] !");
		System.out.println("Finished!");
	}
}
