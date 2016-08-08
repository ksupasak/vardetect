/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class ReconstructSequence {
    
    private ArrayList<MerRead> listMer;
    private ArrayList<Long> listAlgnCode;
    private int beginIdxF,lastIdxF,beginIdxB,lastIdxB;
    private long algnCodeF,algnCodeB;
    private long mask = 268435455;
    private long chrF,chrB, posF, posB;
    private String strandNotF,strandNotB;
    private String resultString;
    
    public ReconstructSequence(){
        this.listMer = new ArrayList();
        this.listAlgnCode = new ArrayList();
    }
    
    public ReconstructSequence(ArrayList<MerRead> inListMer, long inAlgnCodeF, long inAlgnCodeB, int inBeginIdxF, int inLastIdxF, int inBeginIdxB, int inLastIdxB){
        this.listMer = new ArrayList();
        this.listMer.equals(inListMer);
        this.listAlgnCode = new ArrayList();
        this.beginIdxF = inBeginIdxF;
        this.lastIdxF = inLastIdxF;
        this.beginIdxB = inBeginIdxB;
        this.lastIdxB = inLastIdxB;
        this.algnCodeF = inAlgnCodeF;
        this.algnCodeB = inAlgnCodeB;
        
        
                
        this.chrF = ((long)this.algnCodeF)>>29;
        this.posF = ((long)this.algnCodeF&this.mask);
        if(((((long)this.algnCodeF)>>28)&1) == 1){
            this.strandNotF = "+";
        }else if(((((long)this.algnCodeF)>>28)&1) == 0){
            this.strandNotF = "-";
//            this.beginIdxF = inLastIdxF;                        // Strand -  index will switch 
//            this.lastIdxF = inBeginIdxF;
        }
        
        
        
        this.chrB = ((long)this.algnCodeB)>>29;
        this.posB = ((long)this.algnCodeB&this.mask);
        if(((((long)this.algnCodeB)>>28)&1) == 1){
            this.strandNotB = "+";
        }else if(((((long)this.algnCodeB)>>28)&1) == 0){
            this.strandNotB = "-";
//            this.beginIdxB = inLastIdxB;                             // Strand -  index will switch 
//            this.lastIdxB = inBeginIdxB;
        }
        
        System.out.println("Possible pattern : chr"+ chrF + ":chr"+chrB+"\t Strand " + strandNotF + strandNotB + "\tstartIdxF: " + beginIdxF + " StopIdxF: "+lastIdxF+"\tStartIdxB: "+beginIdxB+" StopIdxB: "+lastIdxB);
        this.resultString = "Possible pattern : chr"+ chrF + ":chr"+chrB+"\t Strand " + strandNotF + strandNotB + "\tstartIdxF: " + beginIdxF + " StopIdxF: "+lastIdxF+"\tStartIdxB: "+beginIdxB+" StopIdxB: "+lastIdxB;
        generateReconSequence();
        
    }
    
    public ReconstructSequence(ArrayList<MerRead> inListMer, long inAlgnCodeF, int inBeginIdxF, int inLastIdxF){
        this.listMer = new ArrayList();
        //this.listMer = inListMer;
        this.listMer.equals(inListMer);
        this.listAlgnCode = new ArrayList();
        this.beginIdxF = inBeginIdxF;
        this.lastIdxF = inLastIdxF;
        this.algnCodeF = inAlgnCodeF;      
        
        this.chrF = ((long)this.algnCodeF)>>29;
        this.posF = ((long)this.algnCodeF&this.mask);
        if(((((long)this.algnCodeF)>>28)&1) == 1){
            this.strandNotF = "+";
        }else if(((((long)this.algnCodeF)>>28)&1) == 0){
            this.strandNotF = "-";
//            this.beginIdxF = inLastIdxF;                             // Strand -  index will switch 
//            this.lastIdxF = inBeginIdxF;        
        }
        
        System.out.println("Possible pattern : chr"+ chrF + ":chr"+chrB+"\t Strand " + strandNotF + strandNotB + "\tstartIdxF: " + beginIdxF + " StopIdxF: "+lastIdxF+"\tStartIdxB: "+beginIdxB+" StopIdxB: "+lastIdxB);
        this.resultString = "Possible pattern : chr"+ chrF + ":chr"+chrB+"\t Strand " + strandNotF + strandNotB + "\tstartIdxF: " + beginIdxF + " StopIdxF: "+lastIdxF+"\tStartIdxB: "+beginIdxB+" StopIdxB: "+lastIdxB;
        generateReconSequence();
        
    }
    
    public void addAlignCode(long input){
        this.listAlgnCode.add(input);
    }
    
    public ArrayList getAlignCode(){
        return this.listAlgnCode;
    }
    
    public void addMer(MerRead inMer){
        this.listMer.add(inMer);
    }
    
    public ArrayList getMer(){
        return this.listMer;
    }
    
    public void addListMer(ArrayList<MerRead> inListMer){
        this.listMer.equals(inListMer);
    }
    
    public String getResultString(){
        return this.resultString;
    }
    public void generateReconSequence(){
        for(int i=0;i<listMer.size();i++){
            MerRead dummyMerRead = listMer.get(i);
            int kMer = dummyMerRead.getMeLength();
            long code = dummyMerRead.getMerCode();
            String dnaSequence = SequenceUtil.decodeMer(code, kMer);
            // Combind every mer together to long sequence
            // read index of mer math with start and last 
            // get only first character of string and concatenate all the way to the end
        }
        
    }
}
