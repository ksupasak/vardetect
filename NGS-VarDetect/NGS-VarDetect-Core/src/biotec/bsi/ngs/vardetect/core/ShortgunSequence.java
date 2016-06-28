/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;

/**
 *
 * @author soup
 */
public class ShortgunSequence {
 
    private String seq;
    private String readName;
    ArrayList<MerRead> mers;
    ArrayList<AlignmentData> algns;
    
    public ShortgunSequence(String seq){
        this.seq = seq;
        this.mers = new ArrayList();
        this.algns = new ArrayList();
    }
    
    public void addReadName(String readName){
        this.readName = readName;
    }
    
    public void addMerRead(MerRead mer){
        this.mers.add(mer);
    }
    public void addMerReadByIndex(int idx,MerRead mer){
        this.mers.set(idx, mer);
    }
    
    public int getMerReadSize(){
        return this.mers.size();
    }
    
    public ArrayList<MerRead> getMerRead(){
        return this.mers;
    }
    
    public void addAlignmentData(AlignmentData algn){
        this.algns.add(algn);
    }
    
    public String getReadName(){
        return readName;
    }
    public String getSequence(){
        return seq;
    }
    public int getShortgunLength(){
        return seq.length();
    }
    public void CreateAlignmentData(){
        
    }
}
