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
    private String fullReconSequence;
    private int patternType;                            // Pattern Type indicate the type of posible reconstruct sequence: 0 mean no fusion just sungle align and 1 mean have fusion
    
    public ReconstructSequence(){
        this.listMer = new ArrayList();
        this.listAlgnCode = new ArrayList();
    }
    
    public ReconstructSequence(ArrayList<MerRead> inListMer, long inAlgnCodeF, long inAlgnCodeB, int inBeginIdxF, int inLastIdxF, int inBeginIdxB, int inLastIdxB){
        this.listMer = new ArrayList();
        this.listMer.addAll(inListMer);
        this.listAlgnCode = new ArrayList();
        this.beginIdxF = inBeginIdxF;
        this.lastIdxF = inLastIdxF;
        this.beginIdxB = inBeginIdxB;
        this.lastIdxB = inLastIdxB;
        this.algnCodeF = inAlgnCodeF;
        this.algnCodeB = inAlgnCodeB;
        this.patternType = 1;               // pattern type 1 mean have fusion
        this.fullReconSequence = "";
        
                
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
        this.listMer.addAll(inListMer);
        this.listAlgnCode = new ArrayList();
        this.beginIdxF = inBeginIdxF;
        this.lastIdxF = inLastIdxF;
        this.algnCodeF = inAlgnCodeF;      
        this.patternType = 0;               // pattern type 0 mean no fusion
        this.fullReconSequence = "";
        
        this.chrF = ((long)this.algnCodeF)>>29;
        this.posF = ((long)this.algnCodeF&this.mask);
        if(((((long)this.algnCodeF)>>28)&1) == 1){
            this.strandNotF = "+";
        }else if(((((long)this.algnCodeF)>>28)&1) == 0){
            this.strandNotF = "-";
//            this.beginIdxF = inLastIdxF;                             // Strand -  index will switch (Already switch at ShortgunSequence layer) 
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
    
    public String getFullReconSequence(){
        return this.fullReconSequence;
    }
    
    public void generateReconSequence(){
        System.out.println("Generate Reconstruct Sequence in ReconstructSequence");
        this.fullReconSequence = "";
        
        if(this.patternType == 1){
            for(int i=0;i<this.listMer.size();i++){
                MerRead dummyMerRead = this.listMer.get(i);
                int kMer = dummyMerRead.getMeLength();
                int merIdx = dummyMerRead.getMerIndex();
                long code = (dummyMerRead.getMerCode()>>28); // right shift or do and operation with minus 36 bit
                
//                if((this.beginIdxB-this.lastIdxF) < kMer){
//                    int breakPointF = this.lastIdxF + ((this.beginIdxB - this.lastIdxF)-1);
//                    
//                }else{
//                    int breakPointF = this.lastIdxF + (kMer-1);
//                }   
                
//                int breakPointF = this.lastIdxF + ((this.beginIdxB - this.lastIdxF)-1);
                int breakPointF = this.lastIdxF + (kMer-1);
                int breakPointB = this.lastIdxB + (kMer-1);
                String dnaSequence = SequenceUtil.decodeMer(code, kMer);
//                System.out.println("idx : " +merIdx);
//                if(breakPointF>this.lastIdxB){
//                    if(merIdx >= this.beginIdxF & merIdx <= th){
//                    
//                }else{
                   if(merIdx >= this.beginIdxF & merIdx <= breakPointB){
                    System.out.println("Full Recon Sequence : " + this.fullReconSequence);
                    
                        if (merIdx == breakPointF){
                            this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1)).concat("//F//");
                        }else if(merIdx == this.beginIdxB){
                            this.fullReconSequence = this.fullReconSequence.concat("//B//").concat(dnaSequence.substring(0, 1));
                        }else if(merIdx >= this.lastIdxB && merIdx <= breakPointB){
                            if(merIdx == this.listMer.size()-1){
                                this.fullReconSequence = this.fullReconSequence.concat(dnaSequence).concat("//B//");
                            }else if(merIdx == breakPointB){
                                this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1)).concat("//B//"); 
                            }else{
                                this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
                            }                            
                        }else if(merIdx == this.beginIdxF){
                            this.fullReconSequence = this.fullReconSequence.concat("//F//").concat(dnaSequence.substring(0, 1));
                        }else{
                            this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
                        }
                    
                    } 
//                }
                /*************************/
//                if(merIdx >= this.beginIdxF & merIdx <= this.lastIdxB)
//                    
//                this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
                
                /*if(merIdx >= this.beginIdxF & merIdx <= this.lastIdxB){
//                    System.out.println("Full Recon Sequence : " + this.fullReconSequence);
                    
                    if (merIdx == breakPointF){
                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1)).concat("//F//");
                    }else if(merIdx == this.beginIdxB){
                        this.fullReconSequence = this.fullReconSequence.concat("//B//").concat(dnaSequence.substring(0, 1));
                    }else if(merIdx == this.lastIdxB){
                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence).concat("//B//");
                    }else if(merIdx == this.beginIdxF){
                        this.fullReconSequence = this.fullReconSequence.concat("//F//").concat(dnaSequence.substring(0, 1));
                    }else{
                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
                    }
                    
                }*/
//                else{
//                    if (merIdx == breakPointF){
//                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1)).concat("//~");
//                    }else if (merIdx < breakPointF){
//                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
//                    }
//                }
                
            // Combind every mer together to long sequence
            // read index of mer math with start and last 
            // get only first character of string and concatenate all the way to the end
            }   
        }else if(this.patternType == 0){
             for(int i=0;i<listMer.size();i++){
                MerRead dummyMerRead = listMer.get(i);
                int kMer = dummyMerRead.getMeLength();
                int merIdx = dummyMerRead.getMerIndex();
                long code = (dummyMerRead.getMerCode()>>28);
                int breakPointF = this.lastIdxF + (kMer-1);
                String dnaSequence = SequenceUtil.decodeMer(code, kMer);
//                System.out.println("idx : " +merIdx);
                if(merIdx >= this.beginIdxF & merIdx <= this.lastIdxF){
//                    System.out.println("Full Recon Sequence : " + this.fullReconSequence);
                    if (merIdx == this.lastIdxF){
                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence);
                    }else{
                        this.fullReconSequence = this.fullReconSequence.concat(dnaSequence.substring(0, 1));
                    }
                }
             }
        }
               
        System.out.println("Full Recon Sequence : " + this.fullReconSequence);
        
    }
}
