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
    int readLengthF;

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
    int readLengthB;
    
    char variationType;
    int merLength;
//    int readLength;
    
    long breakPointF;                           // breakpoint value (may be the same as original if it has no change og breakpoint value by other function like recorrectBreakPoint)
    long oriBreakPointF;                        // original breakpoint (this value will not change no effect from any recorectBreakpoint function like recorrectBreakPoint)
    long breakPointB;                           // breakpoint value (may be the same as original if it has no change og breakpoint value by other function like recorrectBreakPoint)
    long oriBreakPointB;                        // original breakpoint (this value will not change no effect from any recorectBreakpoint function like recorrectBreakPoint)
    
    long indelBase;
    String indelType = "";
    
    public Variation(int merLen){
        /**
         * create new object Variation which has merLen and readLen as it has been assign
         * or with default 18 merLen and 100 readLen
         */
        this.merLength = merLen;
//        this.readLength = readLen; 
    }
    
    public Variation(){
        this.merLength = 18;
//        this.readLength = 100; 
    }
    
    public void addFrontPeak(int numChr,long iniPos,long lastPos,int numG,int numY,int numO,int numR,String strand,int iniIdx,String readName,int snpFlag,int iniBackFlag, int inReadLen){
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
        this.readLengthF = inReadLen;
        /**
         * Due to the new alignment protocol, the snpFlag is not working at all so we will comment it out No need to check snpFlag
        */
   
        
//        if(snpFlag == 0){
            if(this.strandF.equals("-")){
                int reverseIniIdx = this.readLengthF-(this.iniIndexF+(this.merLength+this.numMatchF-1));
                this.breakPointF = this.iniPosF+reverseIniIdx;
                this.oriBreakPointF = this.iniPosF+reverseIniIdx;
            }else if(this.strandF.equals("+")){
                this.breakPointF = this.lastPosF+this.iniIndexF;
                this.oriBreakPointF = this.lastPosF+this.iniIndexF;
            }
//        }else{
//            if(this.strandF.equals("-")){
////                int reverseIniIdx = this.readLength-(this.iniBackFlagF+(this.merLength+this.numMatchF-1));
//                this.breakPointF = this.iniPosF+this.iniBackFlagF;   
//            }else if(this.strandF.equals("+")){
//                this.breakPointF = this.iniPosF+this.iniBackFlagF;
//            }
//        }
    }

    public long getOriBreakPointF() {
        return oriBreakPointF;
    }

    public long getOriBreakPointB() {
        return oriBreakPointB;
    }
    
    public void addBackPeak(int numChr,long iniPos,long lastPos,int numG,int numY,int numO,int numR,String strand,int iniIdx,String readName,int snpFlag,int iniBackFlag, int inReadLen){
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
        this.readLengthB = inReadLen;
        
        if(this.strandB != null){
            if(this.strandB.equals("-")){
                int reverseIniIdx = this.readLengthB-(this.iniIndexB+(this.merLength+this.numMatchB-1));
                this.breakPointB = this.lastPosB+reverseIniIdx;
                this.oriBreakPointB = this.lastPosB+reverseIniIdx;
            }else if(this.strandB.equals("+")){
                this.breakPointB = this.iniPosB+this.iniIndexB;
                this.oriBreakPointB = this.iniPosB+this.iniIndexB;
            }
        }
        
        if(this.variationType == 'I'){
            analyzeIndel();
        }else if(this.variationType == 'F'){
            analyzeFusion();
        }
    }
    
    public void addType(char type){
        /**
         * type is represent by one char
         * F = fusion
         * I = indel
         * S = SNP and other un-match between
         * O = others (wasted)
         * T = one Tail
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
    
    public void analyzeIndel(){
        /**
         * This function analyze the indel information
         * 1. indicate type of Indel (insert, delete or large indel)
         * 2. Find indel base in case of insertion and deletion
         * 3. Re-correct break point and number of mer match to relate to real indel base (แก้ไขให้สัมพันธ์กับจำนวน indel base จริงๆ)
         * The re-correction can cause the break point of the pattern not equal to break point from blat (but not too much different)
         * 
         * Implement new way to find indel base in case of strand -+,+- and -- (with strange orientation)
         *  indel base that calculate from alnposF - alnPosB is not correct for above three case
         *  We have to re calculate with another way 
         *          in case -+ ; alnposF will be change to breakPointF - index of breakpointF on original read  (not implement)
         *          in case +- ; alnPosB will be change to breakPointB - index of breakPointB on original read  (not implement)
         *          in case -- ; Both alnPosF and alnPosB will be change to breakPointF - index of breakpointF on original read and breakPointB - index of breakPointB on original read consecutively (Done)
         *          
         *  And also Add case check to distinguish -- with correct order and strange order
         *      we can check by if strand type is -- check the breakPoint => if breakPointF > breakPointB it mean strand -- with correct order
         *      if breakPointF < breakPointB it mean strand -- with incorrect order use special method to find indel base 
         */
        if(this.readNameF.equals("Read12SS00")){
            System.out.println();
        }
        
        this.indelBase = Math.abs(this.iniPosF - this.iniPosB);      // indel base (base on position minus index calculation)
        long breakpointDiff = Math.abs(this.breakPointB - this.breakPointF)-1;
        if(breakpointDiff<=this.readLengthF){
            // It is small Indel (can indicate the number of indel)
            int overallBaseMatch = ((this.numMatchF+this.merLength)-1)+((this.numMatchB+this.merLength)-1);
//            int numBaseF = (this.numMatchF+this.iniIndexF-1);
//            int lastIndexF = numBaseF-1;
//            this.indelBase = this.iniIndexB - lastIndexF;
            
            if(this.strandF.equals("+") && this.strandB.equals("+")){
                /**
                 * for strand ++ 
                 */
                if(this.iniPosB>this.iniPosF){
                    /**
                     * It is deletion
                     * calculate delete base from breakpoint
                     */
                    this.indelType = "delete";
                    int calIndelBase = (int)Math.abs(this.breakPointF - this.breakPointB)-1;
                    
                    if(this.indelBase > calIndelBase){
                        int lastIndexF = (this.iniIndexF+this.numMatchF+this.merLength-2);
                        int overlapBase = Math.abs(this.iniIndexB-lastIndexF)+1;
                        
                        /*
                        Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                        if mer count on front part is too less to Re-correct we will change it to backpart
                        */
                        
                        reCorrectMatchCount(overlapBase);  
                    }

                }else if(this.iniPosB<this.iniPosF){
                    /**
                     * It is insertion
                     * calculate insert base from match index on read
                     */
                    this.indelType = "insert";
                    int numBaseF = (this.numMatchF+this.iniIndexF+this.merLength-1);
                    int lastIndexF = numBaseF-1;
                    int calIndelBase = (this.iniIndexB - lastIndexF)-1;
                    
                    if(this.indelBase > calIndelBase){
                         /*
                        Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                        if mer count on front part is too less to Re-correct we will change it to backpart
                        */
                         
                        if(calIndelBase<0){
                            // case 2 : overlap more than insert base (overlap all insert base and some base on back part)
                            // overlapBase = insert base + number of base that overlap on back part
                            int overlapBase = (int)this.indelBase + Math.abs(calIndelBase);
                            reCorrectMatchCount(overlapBase);
                        }else{
                            // case 1 : overlap on insert base
                            // overlapBase = insertBase - calIndelBase 
                            int overlapBase = (int)this.indelBase - calIndelBase;
                            reCorrectMatchCount(overlapBase);
                        } 
                    }
                }
            }else if(this.strandF.equals("-") && this.strandB.equals("-")){
                /**
                 * for strand -- 
                 * classify criteria change vise versa with strand ++
                 */
                
                if(this.breakPointF < this.breakPointB){
                    // intercept for re calculate indel base if it is -- with wrong order
                    // use new method for calculate indelbase                    
                    
                    int breakPointIndexF = (this.iniIndexF + this.numMatchF + this.merLength -1)-1;      // can use iniIndex directly because we already tranform it to original read index
                    long breakPointAlnPosF = this.breakPointF - breakPointIndexF;                       // calculate position minus index only at break point 
                    
                    int breakPointIndexB = this.iniIndexB;
                    long breakPointAlnPosB = this.breakPointB - breakPointIndexB;
                    
                    // re calculate indelBase with new breakpoint index
                    this.indelBase = Math.abs(breakPointAlnPosF - breakPointAlnPosB);
                }

                if(this.iniPosB>this.iniPosF){
                    /**
                     * It is insertion
                     * calculate insert base from match index on read
                     */
                    this.indelType = "insert";
                    int numBaseF = (this.numMatchF+this.iniIndexF+this.merLength-1);
                    int lastIndexF = numBaseF-1;
                    int calIndelBase = (this.iniIndexB - lastIndexF)-1;

                    if(this.indelBase > calIndelBase){
                        /*
                        Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                        if mer count on front part is too less to Re-correct we will change it to backpart
                        */

                        if(calIndelBase<0){
                            // case 2 : overlap more than insert base (overlap all insert base and some base on back part)
                            // overlapBase = insert base + number of base that overlap on back part
                            int overlapBase = (int)this.indelBase + Math.abs(calIndelBase);
                            reCorrectMatchCount(overlapBase);
                        }else{
                            // case 1 : overlap on insert base
                            // overlapBase = insertBase - calIndelBase 
                            int overlapBase = (int)this.indelBase - calIndelBase;
                            reCorrectMatchCount(overlapBase);
                        } 
                    }
                }else if(this.iniPosB<this.iniPosF){
                    /**
                     * It is deletion
                     * calculate delete base from breakpoint
                     */
//                    this.indelType = "delete";
//                    this.indelBase = Math.abs(this.breakPointF - this.breakPointB)-1;
                    this.indelType = "delete";
                    int calIndelBase = (int)Math.abs(this.breakPointF - this.breakPointB)-1;

                    if(this.indelBase > calIndelBase){
                        int lastIndexF = (this.iniIndexF+this.numMatchF+this.merLength-2);
                        int overlapBase = Math.abs(this.iniIndexB-lastIndexF)+1;

                        /*
                        Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                        if mer count on front part is too less to Re-correct we will change it to backpart
                        */

                        reCorrectMatchCount(overlapBase);  
                    }
                }

            }    
        }else{
            /**
             * It is Large indel
             * We will act with large indel as deletion (No case check for insert or delete like small indel)
             * So, we use deletion method as common.
             */
            this.indelType = "large indel";
            
            if(this.strandF.equals("+") && this.strandB.equals("+")){
                /**
                 * For strand ++ 
                 * Act as deletion
                 * calculate delete base from breakpoint
                 * use indelBase calculated from normal method and have re correct Match Count.
                 */

                int calIndelBase = (int)Math.abs(this.breakPointF - this.breakPointB)-1;

                if(this.indelBase > calIndelBase){
                    int lastIndexF = (this.iniIndexF+this.numMatchF+this.merLength-2);
                    int overlapBase = Math.abs(this.iniIndexB-lastIndexF)+1;

                    /*
                    Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                    if mer count on front part is too less to Re-correct we will change it to backpart
                    */

                    reCorrectMatchCount(overlapBase);  
                }

            }else if(this.strandF.equals("-") && this.strandB.equals("-")){
                /**
                 * For strand -- 
                 * classify criteria change vise versa with strand ++
                 * Act as deletion
                 */
                int calIndelBase = (int)Math.abs(this.breakPointF - this.breakPointB)-1;
                if(this.breakPointF < this.breakPointB){
                    // intercept for re calculate indel base if it is -- with wrong order
                    // Have some evidence that it not suit all the case to re calculate indel base with new method (incase that breakpoint is shift and we use index of shift breakpoint to calculate indelBase. the indelbase that we get is wrong)
                    // So, we set indelBase equal to calIndelBase and no recorrect Match Count.
                    
//                    int breakPointIndexF = (this.iniIndexF + this.numMatchF + this.merLength -1)-1;      // can use iniIndex directly because we already tranform it to original read index
//                    long breakPointAlnPosF = this.breakPointF - breakPointIndexF;                       // calculate position minus index only at break point 
//                    
//                    int breakPointIndexB = this.iniIndexB;
//                    long breakPointAlnPosB = this.breakPointB - breakPointIndexB;
//                    
//                    // re calculate indelBase with new breakpoint index
//                    this.indelBase = Math.abs(breakPointAlnPosF - breakPointAlnPosB);
                    this.indelBase = calIndelBase;
                    
                }else{
                    /**
                     * For strand -- with correct order 
                     * No change with indelBase calculation method and have re correct Match Count.
                     */
                    if(this.indelBase > calIndelBase){
                        int lastIndexF = (this.iniIndexF+this.numMatchF+this.merLength-2);
                        int overlapBase = Math.abs(this.iniIndexB-lastIndexF)+1;

                        /*
                        Re-correct mercount and break point (For front part only because we decide to do on front part if it cannot we do back part instead)
                        if mer count on front part is too less to Re-correct we will change it to backpart
                        */

                        reCorrectMatchCount(overlapBase);  
                    }
                }
            }   
        }
    }

    public void analyzeFusion(){
        /**
         * Find out the the breakpoint is overlap or not
         * Then re-correct breakpoint to get rid off overlap breakPoint 
         */
        
        int lastIdxFront = (this.iniIndexF + this.numMatchF + this.merLength)-2;
        if(lastIdxFront >= this.iniIndexB){
            // Have overlap base around junction. Do reCorrectMatchCount.
            int overlapBase = Math.abs(this.iniIndexB - lastIdxFront)+1;
            reCorrectMatchCount(overlapBase);
        }
        
    }
    
    public long getIndelBase() {
        return indelBase;
    }

    public String getIndelType() {
        return indelType;
    }

    public void setIndelType(String indelType) {
        this.indelType = indelType;
    }
    
    public void reCorrectMatchCount(int overlapBase){
        /**
         * Use for recorrect matchCount relate to indelBase that has been calculate from alnPos (position minus index)
         */
        
        // Actually we can recorrect matchCount on whether front or back part
        // the check case has been use to comfirm that front part has enough number of match count to be re-correct
        // of not it will change to back part instead
        
        
        if(overlapBase < this.numMatchF){ 
            // Re-correct on Front part
            if(this.orangeF > overlapBase){
                this.orangeF = this.orangeF - overlapBase;
            }else if(this.yellowF > overlapBase){
                this.yellowF = this.yellowF - overlapBase;
            }else if(this.greenF > overlapBase){
                this.greenF = this.greenF - overlapBase;
            }else if(this.redF > overlapBase){
                this.redF = this.redF - overlapBase;
            }else{
                /**
                 * Subtract each color by order of orange yellow green red
                 */
                int newOrangeF = this.orangeF - overlapBase;
                if(newOrangeF < 0){
                    this.orangeF = 0;
                    overlapBase = Math.abs(newOrangeF);
                }else{
                    overlapBase = 0;
                }
                
                int newYellowF = this.yellowF - overlapBase;
                if(newYellowF < 0){
                    this.yellowF = 0;
                    overlapBase = Math.abs(newYellowF);
                }else{
                    overlapBase = 0;
                }
                
                int newGreenF = this.greenF - overlapBase;
                if(newGreenF < 0){
                    this.greenF = 0;
                    overlapBase = Math.abs(newGreenF);
                }else{
                    overlapBase = 0;
                }
                
                int newRedF = this.redF - overlapBase;
                if(newRedF < 0){
                    this.redF = 0;
                    overlapBase = Math.abs(newRedF);
                    throw new NumberFormatException();
                }else{
                    overlapBase = 0;
                }    
            }

            this.numMatchF = this.greenF+this.yellowF+this.orangeF+this.redF;
            this.lastPosF = (this.iniPosF + (this.numMatchF+this.merLength-1))-1;
            if(this.strandF.equals("-")){
                int reverseIniIdx = this.readLengthF-(this.iniIndexF+(this.merLength+this.numMatchF-1));
                this.breakPointF = this.iniPosF+reverseIniIdx;   
            }else if(this.strandF.equals("+")){
                this.breakPointF = this.lastPosF+this.iniIndexF;
                System.out.println();
            }   
        }else if(overlapBase < this.numMatchB){
            // Re-correct on Back part
            if(this.orangeB > overlapBase){
                this.orangeB = this.orangeB - overlapBase;
            }else if(this.yellowB > overlapBase){
                this.yellowB = this.yellowB - overlapBase;
            }else if(this.greenB > overlapBase){
                this.greenB = this.greenB - overlapBase;
            }else if(this.redB > overlapBase){
                this.redB = this.redB - overlapBase;
            }else{
                /**
                 * Subtract each color by order of orange yellow green red
                 */
                int newOrangeB = this.orangeB - overlapBase;
                if(newOrangeB < 0){
                    this.orangeB = 0;
                    overlapBase = Math.abs(newOrangeB);
                }else{
                    overlapBase = 0;
                }
                
                int newYellowB = this.yellowB - overlapBase;
                if(newYellowB < 0){
                    this.yellowB = 0;
                    overlapBase = Math.abs(newYellowB);
                }else{
                    overlapBase = 0;
                }
                
                int newGreenB = this.greenB - overlapBase;
                if(newGreenB < 0){
                    this.greenB = 0;
                    overlapBase = Math.abs(newGreenB);
                }else{
                    overlapBase = 0;
                }
                
                int newRedB = this.redB - overlapBase;
                if(newRedB < 0){
                    this.redB = 0;
                    overlapBase = Math.abs(newRedB);
                    throw new NumberFormatException();
                }else{
                    overlapBase = 0;
                }    
            }
            
            this.numMatchB = this.greenB+this.yellowB+this.orangeB+this.redB;
            this.lastPosB = (this.iniPosB + (this.numMatchB+this.merLength-1))-1;
            if(this.strandB != null){
                if(this.strandB.equals("-")){
                    int reverseIniIdx = this.readLengthB-(this.iniIndexB+(this.merLength+this.numMatchB-1));
                    this.breakPointB = this.lastPosB+reverseIniIdx;   
                }else if(this.strandB.equals("+")){
                    this.breakPointB = this.iniPosB+this.iniIndexB;
                }
            }
        }   
    }
    
    public String virtualSequence(){
        /**
         * Create virtual sequence as string
         */
//        if(this.readNameB.equals("HWI-ST840:377:D2GY3ACXX:1:2113:11955:24852/1")&&this.iniPosB == 25556969){
//            System.out.println();
//        }
        int virtualLen = (2*this.readLengthF)+1;    // virtual length has size 2 time of read length plus 1 (plus one is reserve for junction which does not involve with base in virtual length) Ex len = 10 ; virtual = (2*10)+1 = 21 this allow both size have 10 emty slot for base and the middle is junction slot
        int junctionIndex = this.readLengthF;       // index of junction (put sign "|" on this index)
        
        int numBaseMatchF = (this.numMatchF + this.merLength)-1;        // number of base matched (front)
        int numBaseMatchB = (this.numMatchB + this.merLength)-1;        // number of base matched (back)
        int lastIdxF = (this.iniIndexF + numBaseMatchF)-1;
        int idealIniIdxB = lastIdxF + 1;
        int overlapBase = idealIniIdxB - this.iniIndexB;                // if it minus value it mean insertion. Positive is mean overlap and Zero is mean no overlap
        int unMatchBtwJunction = 0;                                     // this variable store the number of un match base that staye in the middle between front and back break point
        
        if(overlapBase < 0){
            // overlap base < 0 is mean insertion or large indel with unmatch . So, no need to change. we set overlapbase to 0.
            unMatchBtwJunction = Math.abs(overlapBase);
            overlapBase = 0;
        }
        numBaseMatchB = numBaseMatchB - overlapBase;                    // calculate new numBaseMatchB if it has overlap (if it not the overlap base is 0. So, no effect at all)
        int newIniIndexB = this.iniIndexB+overlapBase;                  // re calculate iniIndexB if it has overlap (if it not the overlap base is 0. So, no effect at all)
        
//        int virtualIniIndexB = this.iniIndexB*2;                        // the virtual sequence has been create 2 time biger. So for the rigth iniIndexB on the virtual sequence should be multiple by 2 as well
//        int overlapBase = 0;
//        if(virtualIniIndexB < junctionIndex){
//            /**
//             * if back part has overlap with front part. We have to recalculate numBaseMatchB by minus the overlap base out from old numBaseMathB
//             */
//            int lastIdxF = (this.iniIndexF + numBaseMatchF)-1;
//            overlapBase = Math.abs(this.iniIndexB - lastIdxF);
//            numBaseMatchB = numBaseMatchB - overlapBase;
//        }        
        
        int emptySlotPlusUnMatchF = junctionIndex - numBaseMatchF;             // number of empty slot and unMatch slot before match base (front part)
        int numBaseF = numBaseMatchF + this.iniIndexF;
        int emptySlotF = junctionIndex - numBaseF;             // number of empty slot before sequence slot(front part) which include un match base slot and match slot
        
        int remainBaseB = this.readLengthB - (numBaseMatchB + newIniIndexB);      // number of base remain un match on back part;
        int numBaseB = numBaseMatchB + remainBaseB;
        int lastIdxBaseMatchB = junctionIndex + numBaseMatchB;      // last index of match before empty slot (back part)
        int lastIdxBaseB = junctionIndex + numBaseB;
        
        int lastIdxInsertBase = 0;
        int junctionIndexTwo = 0;
        int lastIdxUnMatchBtwJunction = 0;
        if(this.indelType.equals("insert")){
            /**
            * Insertion case => define new variable and adjust old variable for virtualize the insert portion Ex >>>|++++|>>>>  this mean 4 insertion between front and back
            */
           lastIdxInsertBase = junctionIndex + (int)this.indelBase;
           junctionIndexTwo = lastIdxInsertBase + 1;
           lastIdxBaseMatchB = lastIdxBaseMatchB + (int)this.indelBase+1;       // plus 1 to compensate for junction index
           lastIdxBaseB = lastIdxBaseB + (int)this.indelBase+1;                 // plus 1 to compensate for junction index
        }else if(this.indelType.equals("front_oneTail")){
            
        }else if(unMatchBtwJunction!=0){
            /**
             * other case that has un match between junction
             */
            lastIdxUnMatchBtwJunction = junctionIndex + unMatchBtwJunction;
            junctionIndexTwo = lastIdxUnMatchBtwJunction + 1;
            lastIdxBaseMatchB = lastIdxBaseMatchB + unMatchBtwJunction+1;       // plus 1 to compensate for junction index
            lastIdxBaseB = lastIdxBaseB + unMatchBtwJunction+1;                 // plus 1 to compensate for junction index
        }
        
        StringBuilder builder = new StringBuilder();
        
        for(int i=0;i<virtualLen;i++){
            
            if(this.indelType.equals("insert")){
                /**
                 * insert case
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i> junctionIndex && i<=lastIdxInsertBase){
                    builder.append("+");
                }else if(i == junctionIndexTwo){
                    builder.append("|");
                }else if(i>junctionIndexTwo && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }else if(unMatchBtwJunction != 0){
                /**
                 * un match between break point case can happen in large indel or fusion
                 * check with unMatchBtwJunction variable 
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i> junctionIndex && i<=lastIdxUnMatchBtwJunction){
                    builder.append("=");
                }else if(i == junctionIndexTwo){
                    builder.append("|");
                }else if(i>junctionIndexTwo && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }else if(this.indelType.equals("front_oneTail")){
                /**
                 * Case front one Tail 
                 */
                int remainBase = this.readLengthF - numBaseF;
                int lastIndex = junctionIndex+remainBase;
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i>junctionIndex && i<=lastIndex){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
                
            }else if(this.indelType.equals("back_oneTail")){
                /**
                 * Case back one Tail 
                 */
                numBaseMatchB = (this.numMatchB + this.merLength)-1;
                lastIdxBaseMatchB = junctionIndex + numBaseMatchB;
                int remainBaseFront = this.iniIndexB;
                int remainBaseBack = this.readLengthB - (numBaseMatchB+remainBaseFront);
                int firstIndex = junctionIndex-remainBaseFront;
                numBaseB = numBaseMatchB + remainBaseBack;
                lastIdxBaseB = junctionIndex + numBaseB;
                if(i<firstIndex){
                    builder.append(" ");
                }else if(i>=firstIndex && i<junctionIndex){
                    builder.append("=");    
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
                
            }else{
                /**
                 * Normal case + deletion case
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }
//            if(this.indelType.endsWith("insert")){
//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
//                    builder.append("-");
//                }else if(i>lastIdxBaseB){
//                    builder.append(" ");
//                }else{
//                    builder.append("-");
//                }
//                
//            }else{
//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else{
//                    builder.append(" ");
//                }

//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
//                    builder.append("-");
//                }else{
//                    builder.append(" ");
//                }
//            }
            
        }
        
        return builder.toString();
    }
    
    public String oneTailVirtualSequenceFront(){
        /**
         * Create virtual sequence as string 
         * For one tail front event
         */
//        if(this.readNameB.equals("HWI-ST840:377:D2GY3ACXX:1:2113:11955:24852/1")&&this.iniPosB == 25556969){
//            System.out.println();
//        }
        int virtualLen = (2*this.readLengthF)+1;    // virtual length has size 2 time of read length plus 1 (plus one is reserve for junction which does not involve with base in virtual length) Ex len = 10 ; virtual = (2*10)+1 = 21 this allow both size have 10 emty slot for base and the middle is junction slot
        int junctionIndex = this.readLengthF;       // index of junction (put sign "|" on this index)
        
        int numBaseMatchF = (this.numMatchF + this.merLength)-1;        // number of base matched (front)
        int numBaseMatchB = (this.numMatchB + this.merLength)-1;        // number of base matched (back)
        int lastIdxF = (this.iniIndexF + numBaseMatchF)-1;
        int idealIniIdxB = lastIdxF + 1;
        int overlapBase = idealIniIdxB - this.iniIndexB;                // if it minus value it mean insertion. Positive is mean overlap and Zero is mean no overlap
        int unMatchBtwJunction = 0;                                     // this variable store the number of un match base that staye in the middle between front and back break point
        
        if(overlapBase < 0){
            // overlap base < 0 is mean insertion or large indel with unmatch . So, no need to change. we set overlapbase to 0.
            unMatchBtwJunction = Math.abs(overlapBase);
            overlapBase = 0;
        }
        numBaseMatchB = numBaseMatchB - overlapBase;                    // calculate new numBaseMatchB if it has overlap (if it not the overlap base is 0. So, no effect at all)
        int newIniIndexB = this.iniIndexB+overlapBase;                  // re calculate iniIndexB if it has overlap (if it not the overlap base is 0. So, no effect at all)
        
//        int virtualIniIndexB = this.iniIndexB*2;                        // the virtual sequence has been create 2 time biger. So for the rigth iniIndexB on the virtual sequence should be multiple by 2 as well
//        int overlapBase = 0;
//        if(virtualIniIndexB < junctionIndex){
//            /**
//             * if back part has overlap with front part. We have to recalculate numBaseMatchB by minus the overlap base out from old numBaseMathB
//             */
//            int lastIdxF = (this.iniIndexF + numBaseMatchF)-1;
//            overlapBase = Math.abs(this.iniIndexB - lastIdxF);
//            numBaseMatchB = numBaseMatchB - overlapBase;
//        }        
        
        int emptySlotPlusUnMatchF = junctionIndex - numBaseMatchF;             // number of empty slot and unMatch slot before match base (front part)
        int numBaseF = numBaseMatchF + this.iniIndexF;
        int emptySlotF = junctionIndex - numBaseF;             // number of empty slot before sequence slot(front part) which in clude un match base slot and match slot
        
        int remainBaseB = this.readLengthB - (numBaseMatchB + newIniIndexB);      // number of base remain un match on back part;
        int numBaseB = numBaseMatchB + remainBaseB;
        int lastIdxBaseMatchB = junctionIndex + numBaseMatchB;      // last index of match before empty slot (back part)
        int lastIdxBaseB = junctionIndex + numBaseB;
        
        int lastIdxInsertBase = 0;
        int junctionIndexTwo = 0;
        int lastIdxUnMatchBtwJunction = 0;
        if(this.indelType.equals("insert")){
            /**
            * Insertion case => define new variable and adjust old variable for virtualize the insert portion Ex >>>|++++|>>>>  this mean 4 insertion between front and back
            */
           lastIdxInsertBase = junctionIndex + (int)this.indelBase;
           junctionIndexTwo = lastIdxInsertBase + 1;
           lastIdxBaseMatchB = lastIdxBaseMatchB + (int)this.indelBase+1;       // plus 1 to compensate for junction index
           lastIdxBaseB = lastIdxBaseB + (int)this.indelBase+1;                 // plus 1 to compensate for junction index
        }else if(unMatchBtwJunction!=0){
            /**
             * other case that has un match between junction
             */
            lastIdxUnMatchBtwJunction = junctionIndex + unMatchBtwJunction;
            junctionIndexTwo = lastIdxUnMatchBtwJunction + 1;
            lastIdxBaseMatchB = lastIdxBaseMatchB + unMatchBtwJunction+1;       // plus 1 to compensate for junction index
            lastIdxBaseB = lastIdxBaseB + unMatchBtwJunction+1;                 // plus 1 to compensate for junction index
        }
        
        StringBuilder builder = new StringBuilder();
        
        for(int i=0;i<virtualLen;i++){
            
            if(this.indelType.equals("insert")){
                /**
                 * insert case
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i> junctionIndex && i<=lastIdxInsertBase){
                    builder.append("+");
                }else if(i == junctionIndexTwo){
                    builder.append("|");
                }else if(i>junctionIndexTwo && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }else if(unMatchBtwJunction != 0){
                /**
                 * un match between break point case can happen in large indel or fusion
                 * check with unMatchBtwJunction variable 
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i> junctionIndex && i<=lastIdxUnMatchBtwJunction){
                    builder.append("=");
                }else if(i == junctionIndexTwo){
                    builder.append("|");
                }else if(i>junctionIndexTwo && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }else{
                /**
                 * Normal case + deletion case
                 */
                if(i<emptySlotF){
                    builder.append(" ");
                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
                    builder.append("=");    
                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
                    if(this.strandF.equals("+")){
                        builder.append(">");
                    }else if(this.strandF.equals("-")){
                        builder.append("<");
                    }
                }else if(i == junctionIndex){
                    builder.append("|");
                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
                    if(this.strandB.equals("+")){
                        builder.append(">");
                    }else if(this.strandB.equals("-")){
                        builder.append("<");
                    }
                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
                    builder.append("=");
                }else{
                    builder.append(" ");
                }
            }
//            if(this.indelType.endsWith("insert")){
//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
//                    builder.append("-");
//                }else if(i>lastIdxBaseB){
//                    builder.append(" ");
//                }else{
//                    builder.append("-");
//                }
//                
//            }else{
//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else{
//                    builder.append(" ");
//                }

//                if(i<emptySlotF){
//                    builder.append(" ");
//                }else if(i>=emptySlotF && i<emptySlotPlusUnMatchF){
//                    builder.append("-");    
//                }else if(i>=emptySlotPlusUnMatchF && i<junctionIndex){
//                    builder.append("=");
//                }else if(i == junctionIndex){
//                    builder.append("|");
//                }else if(i>junctionIndex && i<=lastIdxBaseMatchB){
//                    builder.append("=");
//                }else if(i>lastIdxBaseMatchB && i<=lastIdxBaseB){
//                    builder.append("-");
//                }else{
//                    builder.append(" ");
//                }
//            }
            
        }
        
        return builder.toString();
    }
    
    public String exportBed12Ref(int extendSize){
        /**
         * this function will calculate the extend right and left position from breakpoint front and back 
         * and export in to bed12 format 
         * 
         * this b3d12 format will be use to cut the DNA sequence from reference to create new refarence that use this variation as a template
         * 
         * EX bed12 format 
         * chr1 0 300 readname 0 + 0 300 0 2 100,100 0,200
         * 
         * chromosome|left most position|right most position|referenceName|score|strand|thickStart|thickStop|itemRGB|num portion that we want to cut(defualt is 2)|size of each portion (extend size)|start pos of each portion
         * Bed format will consider zero based starting position
         * if start = 9 and stop = 20 it will include base at position 10 to 20 (not include 9)
         * 
         * Not finish have to consider stand of variation and bed12 format force us to arrange position in order form low to high
         */
        int leftWingBasePos = 0;
        int rightWingBasePos = 0;
        long leftMostPos = 0;
        long rightMostPos = 0;
        long startLeftWing = 0;
        long startRightWing = 0;
//        if(this.breakPointF > this.breakPointB){
//            leftWingBasePos = this.breakPointB
//        }else{
//            
//        }
        if(this.strandF.equals("+")&&this.strandB.equals("+")){
            leftMostPos = this.breakPointF - extendSize;
            rightMostPos = (this.breakPointB + extendSize)-1;
            startLeftWing = leftMostPos;
            startRightWing = this.breakPointB-1;
        }else if(this.strandF.equals("-")&&this.strandB.equals("-")){
            
        }else if(this.strandF.equals("+")&&this.strandB.equals("-")){
            
        }else if(this.strandF.equals("-")&&this.strandB.equals("+")){
            
        }
        


//        long leftMostPos = this.breakPointF - extendSize;
//        long rightMostPos = (this.breakPointB + extendSize)-1;
//        long startLeftWing = leftMostPos;
//        long startRightWing = this.breakPointB-1;
        String refName = "ref_bf"+this.breakPointF+"_bb"+this.breakPointB;
        int numportion = 2;
        
        String bed12Format = "chr"+this.numChrF+"\t"+leftMostPos+"\t"+rightMostPos+"\t"+refName+"\t0\t+\t"+leftMostPos+"\t"+rightMostPos+"\t0\t"+numportion+"\t"+extendSize+","+extendSize+"\t"+startLeftWing+","+startRightWing;
        return bed12Format;
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
        int iniIndexMatch = this.iniIndexF;                                                 // the correct ini index is index 0
        int lastIndexMatch = (this.iniIndexB + this.numMatchB + this.merLength)-2;          // the correct last index is read length-1  EX it is index 99 if read length is 100
        int unmatchFront = 0;
        int unmatchBack = 0;
        
        /**
         * calculate unmatch front and back
         */
        if(iniIndexMatch>0){
            unmatchFront = iniIndexMatch;
        }
        if(lastIndexMatch<(this.readLengthB-1)){
            unmatchBack = (this.readLengthB-1) - lastIndexMatch;
        }
        /***********************************/
        
        if(this.strandF.equals("+")&&this.strandB.equals("+")){
            iniLeftWing = ((this.iniPosF + this.iniIndexF) - extendSize) - unmatchFront;
            iniRightWing = (this.lastPosB + this.iniIndexB) + 1;
        }else if(this.strandF.equals("-")&&this.strandB.equals("-")){
            int reverseIniIdxF = this.readLengthF-(this.iniIndexF+(this.merLength+this.numMatchF-1));
            int reverseIniIdxB = this.readLengthB-(this.iniIndexB+(this.merLength+this.numMatchB-1));
            iniLeftWing = (this.lastPosF+reverseIniIdxF) + 1;
            iniRightWing = ((this.iniPosB+reverseIniIdxB) - extendSize) - unmatchBack;
        }else if(this.strandF.equals("+")&&this.strandB.equals("-")){
            int reverseIniIdxB = this.readLengthB-(this.iniIndexB+(this.merLength+this.numMatchB-1));
            iniLeftWing = ((this.iniPosF+this.iniIndexF) - extendSize) - unmatchFront;
            iniRightWing = ((this.iniPosB+reverseIniIdxB) - extendSize) - unmatchBack;
        }else if(this.strandF.equals("-")&&this.strandB.equals("+")){
            int reverseIniIdxF = this.readLengthF-(this.iniIndexF+(this.merLength+this.numMatchF-1));
            iniLeftWing = (this.lastPosF+reverseIniIdxF) + 1;
            iniRightWing = (this.lastPosB+this.iniIndexB) + 1;
        }
        
        pointer[0] = iniLeftWing;
        pointer[1] = iniRightWing;
        pointer[2] = iniIndexMatch;
        pointer[3] = lastIndexMatch;
        pointer[4] = unmatchFront;
        pointer[5] = unmatchBack;
        
        return pointer;   
    }
}
