/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
public class RefFaIndex {
    
    private Map<String,ArrayList<Long>> refFaIdx;
    
    public RefFaIndex(String refFaIdxFile) throws FileNotFoundException, IOException{
        RandomAccessFile rbRefIdx = new RandomAccessFile(refFaIdxFile,"r");
        String line;
        String name = "";
        long totalLen = 0;
        long offset = 0;        // offset or in another mean is pointer that point to the number of byte of first base of sequence
        long lineBase = 0;      // number of base on each line (some sequence may cover many of line ex 200 base long may cover 10 line if it wirte 20 base per line)
        long lineWidth = 0;     // number of byte in each line (include the new line) EX 1 byte per base 20 base per line so lineWidth = 21 byte (20 base + newline)
        this.refFaIdx = new LinkedHashMap();     // map that store chromosome name as key and ArrayList of index info as value
        
        while ((line = rbRefIdx.readLine()) != null) {
            String[] linePortion = line.split("\t");
            name = linePortion[0];
            totalLen = Long.parseLong(linePortion[1]);
            offset = Long.parseLong(linePortion[2]);
            lineBase = Long.parseLong(linePortion[3]);
            lineWidth = Long.parseLong(linePortion[4]);
            ArrayList<Long> listIdxInfo = new ArrayList();
            listIdxInfo.add(totalLen);
            listIdxInfo.add(offset);
            listIdxInfo.add(lineBase);
            listIdxInfo.add(lineWidth);
            
            refFaIdx.put(name, listIdxInfo);
        }
    }
    
    public long getTotalLength(String chrName){
        // length of chromosome (number of byte of entire chr)
        
        ArrayList<Long> listIdxInfo = this.refFaIdx.get(chrName);
        
        return listIdxInfo.get(0);
    }
    
    public long getOffSet(String chrName){
        // offset or in another mean is pointer that point to the number of byte of first base of sequence
        
        ArrayList<Long> listIdxInfo = this.refFaIdx.get(chrName);
        
        return listIdxInfo.get(1);
    }
    
    public long getLineBase(String chrName){
        // number of base on each line (some sequence may cover many of line ex 200 base long may cover 10 line if it wirte 20 base per line)
        
        ArrayList<Long> listIdxInfo = this.refFaIdx.get(chrName);
        
        return listIdxInfo.get(2);
    }
    
    public long getLineWidth(String chrName){
        // number of byte in each line (include the new line) EX 1 byte per base 20 base per line so lineWidth = 21 byte (20 base + newline)
        ArrayList<Long> listIdxInfo = this.refFaIdx.get(chrName);
        
        return listIdxInfo.get(3);
    }
    
}
