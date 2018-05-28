/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
public class ReferenceIndex {
    
    private Map<String,Long> refIdxMap;
    private long lineBase;
    
    public ReferenceIndex(){
        refIdxMap = new LinkedHashMap();
    }
    
    public void readfaidx(String filePath) throws FileNotFoundException, IOException{
        
        RandomAccessFile raf = new RandomAccessFile(filePath,"rw");
        String line = null;
        while ((line = raf.readLine()) != null) {
            String[] dummyLine = line.split("\t");
            
            this.refIdxMap.put(dummyLine[0], Long.parseLong(dummyLine[2]));
            this.lineBase = Long.parseLong(dummyLine[3]);
        }
        raf.close();
    }
    
}
