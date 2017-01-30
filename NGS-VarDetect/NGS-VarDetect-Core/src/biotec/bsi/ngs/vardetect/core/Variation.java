/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

/**
 *
 * @author worawich
 */
public class Variation {
    
    int numChrF; 
    long iniPosF; 
    long lastPosF;
    int iniIndexF;
    int greenF;
    int yellowF;
    int orangeF; 
    int redF; 
    int snpFlagF;
    int iniBackFlagF;
    String readNameF;
    String strandF;
    int numMatchF;

    int numChrB; 
    long iniPosB; 
    long lastPosB;
    int iniIndexB;
    int greenB;
    int yellowB;
    int orangeB; 
    int redB; 
    int snpFlagB;
    int iniBackFlagB;
    String readNameB;
    String strandB;
    int numMatchB;
    
    char variationType;
    int merLength;
    int readLength;
    
    long breakPointF;
    long breakPointB;
    
    public Variation(int merLen,int readLen){
        /**
         * create new object Variation which has merLen and readLen as it has been assign
         * or with default 18 merLen and 100 readLen
         */
        this.merLength = merLen;
        this.readLength = readLen; 
    }
    
    public Variation(){
        this.merLength = 18;
        this.readLength = 100; 
    }
    
    public void addFrontPeak(int numChr,long iniPos,long lastPos,int numG,int numY,int numO,int numR,String strand,int iniIdx,String readName,int snpFlag,int iniBackFlag){
        this.numChrF = numChr;
        this.iniPosF = iniPos;
        this.lastPosF=lastPos;
        this.greenF=numG;
        this.yellowF=numY;
        this.orangeF=numO;
        this.redF=numR;
        this.strandF=strand;
        this.iniIndexF=iniIdx;
        this.readNameF=readName;
        this.snpFlagF=snpFlag;
        this.iniBackFlagF=iniBackFlag;
        this.numMatchF = numG+numY+numO+numR;
        
        if(snpFlag == 0){
            if(this.strandF.equals("-")){
                int reverseIniIdx = this.readLength-(this.iniIndexF+(this.merLength+this.numMatchF-1));
                this.breakPointF = this.iniPosF+reverseIniIdx;   
            }else if(this.strandF.equals("+")){
                this.breakPointF = this.lastPosF+this.iniIndexF;
            }
        }else{
            if(this.strandF.equals("-")){
//                int reverseIniIdx = this.readLength-(this.iniBackFlagF+(this.merLength+this.numMatchF-1));
                this.breakPointF = this.iniPosF+this.iniBackFlagF;   
            }else if(this.strandF.equals("+")){
                this.breakPointF = this.iniPosF+this.iniBackFlagF;
            }
        }
    }
    
    public void addBackPeak(int numChr,long iniPos,long lastPos,int numG,int numY,int numO,int numR,String strand,int iniIdx,String readName,int snpFlag,int iniBackFlag){
        this.numChrB = numChr;
        this.iniPosB = iniPos;
        this.lastPosB=lastPos;
        this.greenB=numG;
        this.yellowB=numY;
        this.orangeB=numO;
        this.redB=numR;
        this.strandB=strand;
        this.iniIndexB=iniIdx;
        this.readNameB=readName;
        this.snpFlagB=snpFlag;
        this.iniBackFlagB=iniBackFlag;
        this.numMatchB = numG+numY+numO+numR;
        
       
        if(this.strandB.equals("-")){
            int reverseIniIdx = this.readLength-(this.iniIndexB+(this.merLength+this.numMatchB-1));
            this.breakPointB = this.lastPosB+reverseIniIdx;   
        }else if(this.strandB.equals("+")){
            this.breakPointB = this.iniPosB+this.iniIndexB;
        }
        
    }
    
    public void addType(char type){
        /**
         * type is represent by one char
         * F = fusion
         * I = indel
         * S = SNP and other un-match between
         * O = others (wasted)
         */
        this.variationType = type;
    }
    
    public long calculateBreakPoint(int iniIndex,long alnPos){
        long breakpoint = alnPos-iniIndex;
        return breakpoint;
    }
    
    public long getBreakPointFront(){
        return this.breakPointF;
    }
    
    public long getBreakPointBack(){
        return this.breakPointB;
    }
    
}
