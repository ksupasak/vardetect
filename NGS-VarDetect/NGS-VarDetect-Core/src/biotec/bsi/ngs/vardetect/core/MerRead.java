/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class MerRead {
    
    private long merCode;
    private long merCodeComp;
    private int index,indexComp,merLen;
    private ArrayList<Long> chrPos;             // store align chr and position of this MerRead [chr|position]
    private ArrayList<Long> chrAlgn;            // store align chr and align position(position - index) [chr|align position]
    private ArrayList<Long> chrStrandPos;       // (Arraylist that store all code that align on thisMerRead (can have more than one if it repeat sequence)
    private ArrayList<Long> chrStrandPosComp;   // (Arraylist that store all code that align on thisMerRead compliment strand (can have more than one if it repeat sequence) [chr|Strnd|align Position]
    private ArrayList<Long> chrStrandAlgn;      // (Arraylist that store the position minus mer index of all code that align on thisMerRead (can have more than one if it repeat sequence) [chr|Strnd|align Position]
    private String alignSymbol;
    
    public MerRead(){
        this.chrPos = new ArrayList();
        this.chrAlgn = new ArrayList();
        this.chrStrandPos = new ArrayList();
        this.chrStrandPosComp = new ArrayList();
        this.chrStrandAlgn = new ArrayList();
    
    }
    
    public void addAlignSymbol(String symbol){
        // This function will be call when reconstruct
        
        this.alignSymbol = symbol;      // alignSymbol = 1 mean align | alignSymbol = 0 mean not align
    }
    
    public void addMatchResultStrand(long mer,long[] pos, int idx, long chrNumber,int inMerLen){
        
        this.merLen = inMerLen;
        long[] code = pos;
        this.merCode = mer;
        this.index = idx;
        
//        if(pos!=null&&pos.length>0){
//            int len = pos.length;
//            if(pos[0] > 0){
//                for(int i=0;i<len;i++){
//                    code[i] = (chrNumber<<29)+pos[i]; // leftshift 29 bit because pos is 29 bit length [strand:position]
////                    System.out.println();
////                    System.out.println("This is strand check (must be 1) : " + ((code[i]>>28)&1));
////                    System.out.println();
//                    this.chrStrandPos.add(code[i]);
//                    
////                    System.out.println("***** Do Add strand+ round chr " + chrNumber + "*****" );
//                }
//            }
//        }                
    }
    
    public void addMatchResultStrandCompliment(long mer,long[] pos, int idx, long chrNumber,int inMerLen){
        
        this.merLen = inMerLen;
        long[] code = pos;
        this.merCodeComp = mer;
        this.indexComp = idx;
        
        if(pos!=null&&pos.length>0){
            int len = pos.length;
            if(pos[0] > 0){
                for(int i=0;i<len;i++){
                    code[i] = (chrNumber<<29)+pos[i]; // leftshift 29 bit because pos is 29 bit length [strand:position]
//                    System.out.println();
//                    System.out.println("This is strand check (must be 1) : " + ((code[i]>>28)&1));
//                    System.out.println();
                    this.chrStrandPosComp.add(code[i]);
//                    System.out.println("***** Do Add strand- round chr " + chrNumber + "*****" );
                }
            }
        }                
    }
    
    // pos must be a compose number of chr:position and array of long (or long[])
    public void addMatchResult(long mer,long[] pos, int idx, long chrNumber){       // Currently not use
        
        long[] code = pos;
        this.merCode = mer;
        this.index = idx;
        
        if(pos!=null&&pos.length>0){
            int len = pos.length;
            if(pos[0] > 0){
                for(int i=0;i<len;i++){
                    code[i] = (chrNumber<<28)+pos[i]; 
                    this.chrPos.add(code[i]);
                }
            }
        }                
    }
    
    public void createAlignmentResultStrand(){                                      // Currently not use
        String strandNot = "no value";
        for(int i =0;i<this.chrStrandPos.size();i++){
            
            long newAlignResult = (chrStrandPos.get(i)-this.index);
            if(((chrStrandPos.get(i)>>28)&1) == 1){
                strandNot = "+";
            }else if(((chrStrandPos.get(i)>>28)&1) == 0){
                strandNot = "-";
            }        
//            System.out.println();
//            System.out.println("Check chrPos " + chrStrandPos.get(i));
//            System.out.println("Check strand : " + strandNot);
//            System.out.println("Check chrPos tranform to chrnumber" + (chrStrandPos.get(i)>>29));
//            System.out.println();
            chrStrandAlgn.add(newAlignResult);
            
        }        
    }
    
    public void createAlignmentResultStrandV2(){ // Consider both compliment and not compliment
        
        /* First round non compliment */
        String strandNot = "no value";
        for(int i =0;i<this.chrStrandPos.size();i++){
            
            long newAlignResult = (chrStrandPos.get(i)-this.index);     // alignment positiom code come from chrStrandPos - index (strand +) [different strand have their own aligncode and index]
            if(((chrStrandPos.get(i)>>28)&1) == 1){
                strandNot = "+";
            }else if(((chrStrandPos.get(i)>>28)&1) == 0){
                strandNot = "-";
            }        
//            System.out.println();
//            System.out.println("Check chrPos " + chrStrandPos.get(i));
//            System.out.println("Check strand : " + strandNot);
//            System.out.println("Check chrPos tranform to chrnumber" + (chrStrandPos.get(i)>>29));
//            System.out.println();
            chrStrandAlgn.add(newAlignResult);
            
        }
        
        /* Second round compliment */
        strandNot = "no value";
        for(int i =0;i<this.chrStrandPosComp.size();i++){
            
            long newAlignResult = (chrStrandPosComp.get(i)-this.indexComp); // alignment positiom code come from chrStrandPosComp - indexComp (strand -) [different strand have their own aligncode and index]
            if(((chrStrandPosComp.get(i)>>28)&1) == 1){
                strandNot = "+";
            }else if(((chrStrandPosComp.get(i)>>28)&1) == 0){
                strandNot = "-";
            }        
//            System.out.println();
//            System.out.println("Check chrPos " + chrStrandPos.get(i));
//            System.out.println("Check strand : " + strandNot);
//            System.out.println("Check chrPos tranform to chrnumber" + (chrStrandPos.get(i)>>29));
//            System.out.println();
            chrStrandAlgn.add(newAlignResult);
            
        } 
    }
    
    public void createAlignmentResult(){                        // Currently not use
        for(int i =0;i<this.chrPos.size();i++){
            
            long newAlignResult = (chrPos.get(i)-this.index);
            System.out.println();
            System.out.println("Check chrPos " + chrPos.get(i));
            System.out.println("Check chrPos transform to chrnumber" + (chrPos.get(i)>>28));
            System.out.println();
            chrAlgn.add(newAlignResult);
            
        }        
    }
    
    public ArrayList<Long> getAlignmentResultStrand(){
        return this.chrStrandAlgn;
    }
    
    public ArrayList<Long> getAlignmentResult(){                // Currently not use
        return this.chrAlgn;
    }

    public long getMerCode(){
        return this.merCode;
    }
    
    public long getMerCodeComp(){
        return this.merCodeComp;
    }
    
    public int getMerIndex(){
        return this.index;
    }
    
    public int getMerCompIndex(){
        return this.indexComp;
    }
    
    public int getMeLength(){
        return this.merLen;
    }
}
