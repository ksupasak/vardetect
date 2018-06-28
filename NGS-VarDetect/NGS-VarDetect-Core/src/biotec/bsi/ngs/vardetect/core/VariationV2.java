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
    int iniIndexF;
    int numBaseMatchF;
    int lastIndexF;
    int iniPosF;          // first base real position on read of front pattern (if read is minus strand first base will have higher position than break point)
    int lastPosF;           //last base real position   (if minus strand lastPosF > iniPosF)
    int iniIndexB;
    int numBaseMatchB;
    int lastIndexB;
    int iniPosB;          // first base real position
    int lastPosB;           // last base real position (if minus strand lastposB < iniPosB)
    byte flagUnmap;
    
    public VariationV2(){
        this.numBaseMatchF = 0;
        this.numBaseMatchB = 0;
    }

    public byte getFlagUnmap() {
        return flagUnmap;
    }

    public void setFlagUnmap(byte flagUnmap) {
        this.flagUnmap = flagUnmap;
    }

    public int getReadID() {
        return readID;
    }

    public String getReadSeq() {
        return readSeq;
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
        this.breakpointB = breakpointB;   // add +1 for create new referene purpose which cut the rightwing part uncorrect the cause of this mistake still unknown        
        
        // After set last require information break point back it will generate 
        // iniPos Front and lastPos back
        // iniIndexF and iniIndexB
        
//        if(this.readID==10666){
//            System.out.println();
//        }

        /**
         * Recorrect all breakpoint information
         * (temporary fixed)
         */
        this.breakpointIndexB = this.breakpointIndexB+1;
        int breakpointIndexFInv=0;
        int breakpointIndexBInv=0;
        if(this.strandF == 0 &&this.strandB == 0){
            this.breakpointB = this.alignPosB+this.breakpointIndexB;
        }else if(this.strandF == 1 &&this.strandB==1){
            breakpointIndexFInv = (this.readLen - this.breakpointIndexF)-1;
            breakpointIndexBInv = (this.readLen - this.breakpointIndexB)-1;
            this.breakpointF = this.alignPosF + breakpointIndexFInv;
            this.breakpointB = this.alignPosB + breakpointIndexBInv;
        }else if(this.strandF==0&&this.strandB==1){
            breakpointIndexBInv = (this.readLen - this.breakpointIndexB)-1;
            this.breakpointB = this.alignPosB + breakpointIndexBInv;
        }else if(this.strandF==1&&this.strandB==0){
            breakpointIndexFInv = (this.readLen - this.breakpointIndexF)-1;
            this.breakpointF = this.alignPosF + breakpointIndexFInv;
        }
        /*************************************/

        
        boolean flagA = false;
        int sCount = 0; // number of s character in merCollection 
        /**
         * loop find firstA
         */
        boolean fullFlagF = true;
        for(int i=this.breakpointIndexF;i>=0;i--){
            if(this.merCollection.charAt(i) == 'A' ||this.merCollection.charAt(i) == 'a'){
                this.iniIndexF = i;
            }
            
//            if(this.merCollection.charAt(i) == 'A' ||this.merCollection.charAt(i) == 'a'){
//                this.numBaseMatchF++; 
//                sCount=0;
//            }else{
//                sCount++;
//                if(sCount>4){
//                    this.iniIndexF = i+5;
//                    fullFlagF=false;
//                    break;
//                }        
////                this.iniIndexF = i+1;
////                fullFlagF=false;
////                break;
//            }
        }
//        if(fullFlagF==true){
//            this.iniIndexF = 0;
//        }
        
        /**
         * loop find lastB
         */
        sCount=0;   // reset Count
        boolean fullFlagB = true;
        for(int i=this.breakpointIndexB;i<this.merCollection.length();i++){
            if(this.merCollection.charAt(i) == 'B' || this.merCollection.charAt(i) == 'b'){
                this.lastIndexB = i;
            }

//            if(this.merCollection.charAt(i) == 'B' || this.merCollection.charAt(i) == 'b'){
//                this.numBaseMatchB++;
//                sCount=0;
//            }else{
//                sCount++;
//                if(sCount>4){
//                    this.lastIndexB = i-5;
//                    fullFlagF=false;
//                    break;
//                }
////                this.lastIndexB = i-1;
////                fullFlagB=false;
////                break;
//            }
        }
//        if(fullFlagB==true){
//            this.lastIndexB = this.merCollection.length()-1;
//        }
        
//        for(int i=0;i<this.merCollection.length();i++){
//            if(this.merCollection.charAt(i) == 'A' ||this.merCollection.charAt(i) == 'a'){
//                
//                if(flagA == false){
//                    this.iniIndexF = i;
//                    flagA = true;
//                }
//                this.numBaseMatchF++;    
//            }else if(this.merCollection.charAt(i) == 'B' || this.merCollection.charAt(i) == 'b'){
//
//                this.lastIndexB = i;
//                this.numBaseMatchB++;
//            }
//        }
        this.lastIndexF = this.breakpointIndexF;
        this.iniIndexB = this.breakpointIndexB;
        
        this.numBaseMatchF = (this.lastIndexF - this.iniIndexF) + 1;
        this.numBaseMatchB = (this.lastIndexB - this.iniIndexB) + 1;
        
        if(this.strandF == 0 &&this.strandB == 0){
            this.iniPosF = this.breakpointF-this.numBaseMatchF+1;
            this.lastPosF = this.breakpointF;
            this.iniPosB = this.breakpointB;
            this.lastPosB = this.breakpointB+this.numBaseMatchB-1;
        }else if(this.strandF == 1 &&this.strandB==1){
            this.iniPosF = this.breakpointF+this.numBaseMatchF-1;
            this.lastPosF = this.breakpointF;
            this.iniPosB = this.breakpointB;
            this.lastPosB = this.breakpointB-this.numBaseMatchB+1;
        }else if(this.strandF==0&&this.strandB==1){
            this.iniPosF = this.breakpointF-this.numBaseMatchF+1;
            this.lastPosF = this.breakpointF;
            this.iniPosB = this.breakpointB;
            this.lastPosB = this.breakpointB-this.numBaseMatchB+1;
        }else if(this.strandF==1&&this.strandB==0){
            this.iniPosF = this.breakpointF+this.numBaseMatchF-1;
            this.lastPosF = this.breakpointF;
            this.iniPosB = this.breakpointB;
            this.lastPosB = this.breakpointB+this.numBaseMatchB-1;
        }
    }

    public int getBreakpointF() {
        return breakpointF;
    }

    public int getBreakpointB() {
        return breakpointB;
    }

    public int getBreakpointIndexF() {
        return breakpointIndexF;
    }

    public int getBreakpointIndexB() {
        return breakpointIndexB;
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
    
    @Override       // overide equals for check smilarity of two variationV2 object for implement remove duplication read
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
    
    public long[] getSignatureForCreateRef(int extendSize){
        /**
         * This function will create key signature value for create reference
         * create initial pointer of left and right wing
         * return in long[] first element is left wing pointer  second element is right wing pointer  third is iniIndexMatch(first index that match on read) and fourth is lastIndexMatch(last index that match on read)
         * fifth is unmatchFront and sixth is unmatchBack
         */
        long[] pointer = new long[6];           // array of long to store iniLeftWing iniRightWing iniIndexMatch(first index that match on read) and lastIndexMatch(last index that match on read)
        long iniLeftWing = 0;
        long iniRightWing = 0;
//        if(this.breakPointF > this.breakPointB){
//            leftWingBasePos = this.breakPointB
//        }else{
//            
//        }
        int iniIndexMatch = this.iniIndexF;                                                 // index start at 0 so if read has 150 long that last index will be 149
        int lastIndexMatch = this.lastIndexB;           
        int unmatchFront = 0;
        int unmatchBack = 0;
        
        /**
         * calculate unmatch front and back
         */
        if(this.iniIndexF > 0){
            unmatchFront = iniIndexMatch;
        }
        if(lastIndexMatch<(this.readLen-1)){
            unmatchBack = ((this.readLen-1) - lastIndexMatch)-1;
        }
        /***********************************/
        
        if(this.strandF == 0 &&this.strandB == 0){
            iniLeftWing = (this.iniPosF - extendSize) - unmatchFront;
            iniRightWing = this.lastPosB + 1;
        }else if(this.strandF == 1 &&this.strandB==1){
            iniLeftWing = this.iniPosF + 1;
            iniRightWing = (this.lastPosB - extendSize) - unmatchBack;
        }else if(this.strandF==0&&this.strandB==1){            
            iniLeftWing = (this.iniPosF - extendSize) - unmatchFront;
            iniRightWing = (this.lastPosB - extendSize) - unmatchBack;
        }else if(this.strandF==1&&this.strandB==0){            
            iniLeftWing = this.iniPosF + 1;
            iniRightWing = this.lastPosB + 1;
        }
        
        pointer[0] = iniLeftWing;
        pointer[1] = iniRightWing;
        pointer[2] = iniIndexMatch;
        pointer[3] = lastIndexMatch;
        pointer[4] = unmatchFront;
        pointer[5] = unmatchBack;
        
        return pointer;   
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
    
//    @Override
//    public int compareTo(VariationV2 inVar) {
//        int compareCov = ((SVGroup)compareSVGroup).getNumCoverage();
//        /* For Ascending order*/
////        return this.studentage-compareage;
//        /* For Descending order do like this */
//        return compareCov-getNumCoverage();
//    }
}
