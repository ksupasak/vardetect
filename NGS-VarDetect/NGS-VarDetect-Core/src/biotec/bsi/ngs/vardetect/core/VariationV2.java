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
public class VariationV2 {
    
    int readID;
    int merCoverage;
    byte orientationCode;
    int numOverlapMer;
    int numMerF;
    int numMerB;
    long posCodeF;      // raw position (can be translate to chromosome and position on reference)
    long posCodeB;      // raw position (can be translate to chromosome and position on reference) 
    byte strandF;   // 0 = + and 1 = -
    byte strandB;
    String chrF;
    String chrB;
    int alignPosF;      // A front part offset position (minus position) 
    int alignPosB;      // A back part offset position (minus position) 
    String readSeq;
    String refF;
    String refB;
    String merProfile;
    String merCollection;
    int breakpointIndexF;       // index on read of last base of front pattern
    int breakpointIndexB;       // index on read of base before first base of back pattern (Must plus one to get index of first base of back pattern)  
    int breakpointF;
    int breakpointB;
    String middleBase;
    int readLen;
    int numBaseMatchF;
    int numBaseMatchB;
    
    
            
    
    public VariationV2(){
        
    }

    public void setReadID(int readID) {
        this.readID = readID;
    }

    public void setMerCoverage(int merCoverage) {
        this.merCoverage = merCoverage;
    }

    public void setOrientationCode(byte orientationCode) {
        this.orientationCode = orientationCode;
    }

    public void setNumOverlapMer(int numOverlapMer) {
        this.numOverlapMer = numOverlapMer;
    }

    public void setNumMerF(int numMerF) {
        this.numMerF = numMerF;
    }

    public void setNumMerB(int numMerB) {
        this.numMerB = numMerB;
    }

    public void setPosCodeF(long posCodeF) {
        this.posCodeF = posCodeF;
    }

    public void setPosCodeB(long posCodeB) {
        this.posCodeB = posCodeB;
    }

    public void setStrandF(byte strandF) {
        this.strandF = strandF;
    }

    public void setStrandB(byte strandB) {
        this.strandB = strandB;
    }

    public void setChrF(String chrF) {
        this.chrF = chrF;
    }

    public void setChrB(String chrB) {
        this.chrB = chrB;
    }

    public void setAlignPosF(int alignPosF) {
        this.alignPosF = alignPosF;
    }

    public void setAlignPosB(int alignPosB) {
        this.alignPosB = alignPosB;
    }

    public void setReadSeq(String readSeq) {
        this.readSeq = readSeq;
        this.readLen = readSeq.length();
    }

    public void setRefF(String refF) {
        this.refF = refF;
    }

    public void setRefB(String refB) {
        this.refB = refB;
    }

    public void setMerProfile(String merProfile) {
        this.merProfile = merProfile;
    }

    public void setMerCollection(String merCollection) {
        this.merCollection = merCollection;
    }

    public void setBreakpointIndexF(int breakpointIndexF) {
        this.breakpointIndexF = breakpointIndexF;
//        this.numBaseMatchF = breakpointIndexF;
    }

    public void setBreakpointIndexB(int breakpointIndexB) {
        this.breakpointIndexB = breakpointIndexB;
    }

    public void setBreakpointF(int breakpointF) {
        this.breakpointF = breakpointF;
    }

    public void setBreakpointB(int breakpointB) {
        this.breakpointB = breakpointB;
    }

    public int getBreakpointF() {
        return breakpointF;
    }

    public int getBreakpointB() {
        return breakpointB;
    }

    public byte getStrandF() {
        return strandF;
    }

    public byte getStrandB() {
        return strandB;
    }

    public int getAlignPosF() {
        return alignPosF;
    }

    public int getAlignPosB() {
        return alignPosB;
    }

    public String getChrF() {
        return chrF;
    }

    public String getChrB() {
        return chrB;
    }

    public String getMiddleBase() {
        return middleBase;
    }

    public void setMiddleBase(String middleBase) {
        this.middleBase = middleBase;
    }

    public long getPosCodeF() {
        return posCodeF;
    }

    public long getPosCodeB() {
        return posCodeB;
    }
    
    @Override
    public boolean equals(Object obj){
        boolean isEqual = false;
        VariationV2 anotherObj = (VariationV2)obj;
        if(obj != null && obj instanceof VariationV2){
            
            if(this.alignPosF == anotherObj.getAlignPosF() && this.alignPosB == anotherObj.getAlignPosB() && this.strandF == anotherObj.strandF && this.strandB == anotherObj.strandB) {
                isEqual=true;
            }  
        }
        return isEqual;
    }
    
    @Override
    public int hashCode(){
        return this.alignPosF;
    }
    
    @Override
    public String toString(){
        return readID+"\t"+merCoverage+"\t"+orientationCode+"\t"+numOverlapMer+"\t"+numMerF+"\t"+numMerB+"\t"+posCodeF+":"+strandF+"\t"+posCodeB+":"+strandB+"\t"
                +chrF+":"+alignPosF+"\t"+chrB+":"+alignPosB+"\t"+readSeq+"\t"+refF+"\t"+refB+"\t"+merProfile+"\t"+merCollection+"\t"+breakpointIndexF+"\t"+breakpointIndexB+"\t"
                +chrF+":"+breakpointF+"\t"+chrB+":"+breakpointB+"\t"+middleBase;
    }
}
