/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
//public class SVGroup implements Comparable<SVGroup> {
public class SVGroup {
    /**
     * Object that use to store group of object variationV2 which has the same structure variant pattern
     */
    private ArrayList<VariationV2> varList;
    private Map<String,Long> refIndex;
    private String svType;
    private int RPB;        // reference break point back
    private int RPF;        // reference break point front
    private int APB;
    private int APF;
    private long rawPosF;
    private long rawPosB;
    private String chrF;
    private String chrB;
    private long numChrF;
    private long numChrB;
    private byte strandF;
    private byte strandB;
    private int numCoverage;
    private byte svTypeCode;    // 0 is tanden, 1 is indel, 2 is intraTrans, 3 is interTrans, 4 is unclassify    
    private boolean ppFlag;     // ++ strad flag
    private boolean mmFlag;     // -- strand flag
    private boolean pmFlag;     // +- strand flag
    private boolean mpFlag;     // -+ strand flag
    private byte numStrandPattern;      // 1 mean this group has one strand pattern , 2 mean has two Strand Pattern (Ex ++ --) ans so on
    private boolean identityFlag;
    private long frontCode;         // [numchr:RPF]
    private long backCode;
    private Annotation annoF;       // object store annotation information of front part
    private Annotation annoB;       // object store annotation information of back part
    private boolean sameChrFlag;
    private boolean diffChrFlag;
    private ArrayList<Double> sameChrEnclideanDistance;
    private ArrayList<Double> diffChrEuclideanDistance;
    private boolean sameChrEnclideanDistanceCompleteFlag;
    private boolean diffChrEuclideanDistanceCompleteFlag;
    private long euCalFrontCode;
    private long euCalBackCode;
    private boolean tandemFlag;
    private boolean deleteFlag;
    private boolean intraTransFlag;
    private boolean interTransFlag;
    private boolean chimericFlag;
    private boolean alienFlag;
    
    public SVGroup(){
        varList = new ArrayList();        
        this.ppFlag=false;
        this.mmFlag=false;
        this.pmFlag=false;
        this.mpFlag=false;
        this.numStrandPattern=0;
        this.refIndex=new LinkedHashMap();
        this.identityFlag=false;
        this.annoF = new Annotation();      
        this.annoB = new Annotation();
        this.sameChrEnclideanDistance = new ArrayList();
        this.diffChrEuclideanDistance = new ArrayList();
        this.sameChrEnclideanDistanceCompleteFlag = false;
        this.diffChrEuclideanDistanceCompleteFlag = false;
        this.tandemFlag = false;
        this.deleteFlag = false;
        this.intraTransFlag = false;
        this.interTransFlag = false;
        this.chimericFlag = false;
        this.alienFlag = false;
    }
    
    public void addVariationV2(VariationV2 inVar){
        if(!this.varList.contains(inVar)){
            this.varList.add(inVar);
            
            if(inVar.getStrandF()==0 && inVar.getStrandB()==0 && this.ppFlag == false){
                this.ppFlag = true;
                this.numStrandPattern++;
            }else if(inVar.getStrandF()==1 && inVar.getStrandB()==1 && this.mmFlag==false){
                this.mmFlag = true;
                this.numStrandPattern++;
            }else if(inVar.getStrandF()==0 && inVar.getStrandB()==1 && this.pmFlag==false){
                this.pmFlag = true;
                this.numStrandPattern++;
            }else if(inVar.getStrandF()==1 && inVar.getStrandB()==0 && this.mpFlag==false){
                this.mpFlag = true;
                this.numStrandPattern++;
            }
        } 
    }
    
    public void addRefIndex(Map<String,Long> inRefIndex){
        this.refIndex = inRefIndex;
    }

    public byte getNumStrandPattern() {
        return numStrandPattern;
    }

    public String getSVType(){
//        VariationV2 dummyVar = this.varList.get(0);
//        this.RPF = dummyVar.getBreakpointF();
//        this.RPB = dummyVar.getBreakpointB();
//        this.APF = dummyVar.getAlignPosF();
//        this.APB = dummyVar.getAlignPosB();
//        this.chrF = dummyVar.getChrF();
//        this.chrB = dummyVar.getChrB();
//        this.strandF = dummyVar.getStrandF();
//        this.strandB = dummyVar.getStrandB();
//        this.rawPosF = dummyVar.getPosCodeF();
//        this.rawPosB = dummyVar.getPosCodeB();
        defineIdentity();
        /**
         * This function will classify rough SV type for this SVGroup
         */
        if(this.chrF.equals(this.chrB)){
            /**
             * Same chromosome
             */
            if(this.strandF==0 && this.strandB==0){
                // Strand ++
                if(this.RPB < this.RPF && this.APB < this.APF){
                    this.svType = "tandem";
                    this.svTypeCode=0;
                }else if(this.RPB > this.RPF && this.APB > this.APF){
                    this.svType = "indel";
                    this.svTypeCode=1;
                }else{
                    this.svType = "unclassify";
                    this.svTypeCode=4;
                }
            }else if(this.strandF==1 && this.strandB==1){
                // Strand --
                if(this.RPB > this.RPF && this.APB > this.APF){
                    this.svType = "tandem";
                    this.svTypeCode=0;
                }else if(this.RPB < this.RPF && this.APB < this.APF){
                    this.svType = "indel";
                    this.svTypeCode=1;
                }else{
                    this.svType = "unclassify";
                    this.svTypeCode=4;
                }
            }else if(this.strandF==0 && this.strandB==1){
                // Strand +-
                if(this.RPB < this.RPF && this.APB < this.APF){
                    this.svType = "intraTrans";
                    this.svTypeCode=2;
                }else if(this.RPB > this.RPF && this.APB > this.APF){
                    this.svType = "intraTrans";
                    this.svTypeCode=2;
                }else{
                    this.svType = "unclassify";
                    this.svTypeCode=4;
                }
            }else if(this.strandF==1 && this.strandB==0){
                // Strand -+
                if(this.RPB < this.RPF && this.APB < this.APF){
                    this.svType = "intraTrans";
                    this.svTypeCode=2;
                }else if(this.RPB > this.RPF && this.APB > this.APF){
                    this.svType = "intraTrans";
                    this.svTypeCode=2;
                }else{
                    this.svType = "unclassify";
                    this.svTypeCode=4;
                }
            }
            
        }else{
            /**
             * different Chromosome
             */
            this.svType = "interTrans";
            this.svTypeCode=3;
            
        }
        
        return this.svType;
    }
    
    public String getPreciseSVType(){
//        VariationV2 dummyVar = this.varList.get(0);
//        this.RPF = dummyVar.getBreakpointF();
//        this.RPB = dummyVar.getBreakpointB();
//        this.APF = dummyVar.getAlignPosF();
//        this.APB = dummyVar.getAlignPosB();
//        this.chrF = dummyVar.getChrF();
//        this.chrB = dummyVar.getChrB();
//        this.strandF = dummyVar.getStrandF();
//        this.strandB = dummyVar.getStrandB();
//        this.rawPosF = dummyVar.getPosCodeF();
//        this.rawPosB = dummyVar.getPosCodeB();

        /**
         * This function should be call after passing defineIdentity
         * And complete calculate and add euclidean distance to arraylist
         * 
         * Because it require identity information  of this SVGroup and ArrayList of euclidean distance of this SVGroup with other SVGroup
         * 
         * To precisely classify SV type
         * 
         * not finish
        */
        
        
        /**
         * This function will classify rough SV type for this SVGroup
         */
//        if(this.chrF.equals(this.chrB)){
//            /**
//             * Same chromosome
//             */
//            if(this.strandF==0 && this.strandB==0){
//                // Strand ++
//                if(this.RPB < this.RPF && this.APB < this.APF){
//                    this.svType = "tandem";
//                    this.svTypeCode=0;
//                }else if(this.RPB > this.RPF && this.APB > this.APF){
//                    this.svType = "indel";
//                    this.svTypeCode=1;
//                }else{
//                    this.svType = "unclassify";
//                    this.svTypeCode=4;
//                }
//            }else if(this.strandF==1 && this.strandB==1){
//                // Strand --
//                if(this.RPB > this.RPF && this.APB > this.APF){
//                    this.svType = "tandem";
//                    this.svTypeCode=0;
//                }else if(this.RPB < this.RPF && this.APB < this.APF){
//                    this.svType = "indel";
//                    this.svTypeCode=1;
//                }else{
//                    this.svType = "unclassify";
//                    this.svTypeCode=4;
//                }
//            }else if(this.strandF==0 && this.strandB==1){
//                // Strand +-
//                if(this.RPB < this.RPF && this.APB < this.APF){
//                    this.svType = "intraTrans";
//                    this.svTypeCode=2;
//                }else if(this.RPB > this.RPF && this.APB > this.APF){
//                    this.svType = "intraTrans";
//                    this.svTypeCode=2;
//                }else{
//                    this.svType = "unclassify";
//                    this.svTypeCode=4;
//                }
//            }else if(this.strandF==1 && this.strandB==0){
//                // Strand -+
//                if(this.RPB < this.RPF && this.APB < this.APF){
//                    this.svType = "intraTrans";
//                    this.svTypeCode=2;
//                }else if(this.RPB > this.RPF && this.APB > this.APF){
//                    this.svType = "intraTrans";
//                    this.svTypeCode=2;
//                }else{
//                    this.svType = "unclassify";
//                    this.svTypeCode=4;
//                }
//            }
//            
//        }else{
//            /**
//             * different Chromosome
//             */
//            this.svType = "interTrans";
//            this.svTypeCode=3;
//            
//        }
        
        return this.svType;
    }
    
    public void defineIdentity(){
        /**
         * This function will define identity for this SVGroup
         */
        this.identityFlag = true; // true mean this SV group already define identity
        if(this.numStrandPattern>1){
            VariationV2 selectVar = this.varList.get(0);
            this.strandF = selectVar.getStrandF();
            this.strandB = selectVar.getStrandB();
            
            if(this.strandF == this.strandB){
                /**
                 * this group is ++ and -- pattern
                 * loop to find ++ strand pattern to use as the main identity of this group
                 */
                for(int i=0;i<this.varList.size();i++){
                    VariationV2 dummyVar = this.varList.get(i);
                    this.strandF = dummyVar.getStrandF();
                    this.strandB = dummyVar.getStrandB();
                    if(this.strandF == 0 && this.strandB == 0){
                        this.RPF = dummyVar.getBreakpointF();
                        this.RPB = dummyVar.getBreakpointB();
                        this.APF = dummyVar.getAlignPosF();
                        this.APB = dummyVar.getAlignPosB();
                        this.chrF = dummyVar.getChrF();
                        this.chrB = dummyVar.getChrB();
                        this.strandF = dummyVar.getStrandF();
                        this.strandB = dummyVar.getStrandB();
                        this.rawPosF = dummyVar.getPosCodeF();
                        this.rawPosB = dummyVar.getPosCodeB();
                        this.numChrF = this.refIndex.get(this.chrF);
                        this.numChrB = this.refIndex.get(this.chrB);
                        this.frontCode = (this.numChrF<<32)+this.RPF;
                        this.backCode = (this.numChrB<<32)+this.RPB;
                        this.euCalFrontCode = (this.strandF<<63)+((this.numChrF<<32)+this.RPF);
                        this.euCalBackCode = (this.strandB<<63)+((this.numChrB<<32)+this.RPB);
                        
                        if(numChrF==numChrB){
                            this.sameChrFlag = true;
                            this.diffChrFlag = false;
                        }else{
                            this.sameChrFlag = false;
                            this.diffChrFlag = true;
                        }
                        
                        break;
                    }  
                }
            }else{
                /**
                 * this. group is +- and -+ pattern
                 * Still confused,so we pick the first as the main identity of this group 
                 */
                VariationV2 dummyVar = this.varList.get(0);
                this.RPF = dummyVar.getBreakpointF();
                this.RPB = dummyVar.getBreakpointB();
                this.APF = dummyVar.getAlignPosF();
                this.APB = dummyVar.getAlignPosB();
                this.chrF = dummyVar.getChrF();
                this.chrB = dummyVar.getChrB();
                this.strandF = dummyVar.getStrandF();
                this.strandB = dummyVar.getStrandB();
                this.rawPosF = dummyVar.getPosCodeF();
                this.rawPosB = dummyVar.getPosCodeB();
                this.numChrF = this.refIndex.get(this.chrF);
                this.numChrB = this.refIndex.get(this.chrB);
                this.frontCode = (this.numChrF<<32)+this.RPF;
                this.backCode = (this.numChrB<<32)+this.RPB;
                this.euCalFrontCode = (this.strandF<<63)+((this.numChrF<<32)+this.RPF);
                this.euCalBackCode = (this.strandB<<63)+((this.numChrB<<32)+this.RPB);
                
                if(numChrF==numChrB){
                    this.sameChrFlag = true;
                    this.diffChrFlag = false;
                }else{
                    this.sameChrFlag = false;
                    this.diffChrFlag = true;
                }
            }
            
            
        }else{
            VariationV2 dummyVar = this.varList.get(0);
            this.RPF = dummyVar.getBreakpointF();
            this.RPB = dummyVar.getBreakpointB();
            this.APF = dummyVar.getAlignPosF();
            this.APB = dummyVar.getAlignPosB();
            this.chrF = dummyVar.getChrF();
            this.chrB = dummyVar.getChrB();
            this.strandF = dummyVar.getStrandF();
            this.strandB = dummyVar.getStrandB();
            this.rawPosF = dummyVar.getPosCodeF();
            this.rawPosB = dummyVar.getPosCodeB();
            this.numChrF = this.refIndex.get(this.chrF);
            this.numChrB = this.refIndex.get(this.chrB);
            this.frontCode = (this.numChrF<<32)+this.RPF;
            this.backCode = (this.numChrB<<32)+this.RPB;
            this.euCalFrontCode = (this.strandF<<63)+((this.numChrF<<32)+this.RPF);
            this.euCalBackCode = (this.strandB<<63)+((this.numChrB<<32)+this.RPB);
            
            if(numChrF==numChrB){
                this.sameChrFlag = true;
                this.diffChrFlag = false;
            }else{
                this.sameChrFlag = false;
                this.diffChrFlag = true;
            }
        }
        
    }
    
    public int getNumCoverage() {
        numCoverage = this.varList.size();
        return numCoverage;
    }
    
    public static Comparator<SVGroup> CoverageComparator = new Comparator<SVGroup>() {

	public int compare(SVGroup s1, SVGroup s2) {
            int s1Coverage = s1.getNumCoverage();
            int s2Coverage = s2.getNumCoverage();
	   //ascending order (low to high)
            //return s1Coverage-s2Coverage;

	   //descending order   (high to low)
            return s2Coverage-s1Coverage;
        }
    };
    
    public static Comparator<SVGroup> FrontBreakPointCodeComparator = new Comparator<SVGroup>() {

	public int compare(SVGroup s1, SVGroup s2) {
            if(s1.identityFlag==false){
                s1.defineIdentity();
            }
            
            if(s2.identityFlag==false){
                s2.defineIdentity();
            }
            
	   //ascending order (low to high)
            if(s1.frontCode < s2.frontCode){
                return -1;
            }else if(s1.frontCode > s2.frontCode){
                return 1;
            }else{
                return 0;
            }
            
	   //descending order (high to low)
//            if(s1.frontCode < s2.frontCode){
//                return 1;
//            }else if(s1.frontCode > s2.frontCode){
//                return -1;
//            }else{
//                return 0;
//            }
        }
    };
    
    public static Comparator<SVGroup> BackBreakPointCodeComparator = new Comparator<SVGroup>() {
        public int compare(SVGroup s1, SVGroup s2) {
            if(s1.identityFlag==false){
                s1.defineIdentity();
            }
            
            if(s2.identityFlag==false){
                s2.defineIdentity();
            }
            
	   //ascending order (low to high)
            if(s1.backCode < s2.backCode){
                return -1;
            }else if(s1.backCode > s2.backCode){
                return 1;
            }else{
                return 0;
            }
            
	   //descending order (high to low)
//            if(s1.frontCode < s2.frontCode){
//                return 1;
//            }else if(s1.frontCode > s2.frontCode){
//                return -1;
//            }else{
//                return 0;
//            }
        }
    };
    
//    @Override
//    public int compareTo(SVGroup compareSVGroup) {
//        int compareCov = ((SVGroup)compareSVGroup).getNumCoverage();
//        /* For Ascending order*/
////        return this.studentage-compareage;
//        /* For Descending order do like this */
//        return compareCov-getNumCoverage();
//    }

    public boolean isIdentityFlag() {
        return identityFlag;
    }

    public Map<String, Long> getRefIndex() {
        return refIndex;
    }

    public void setSvType(String svType) {
        this.svType = svType;
    }

    public String getSvType() {
        return svType;
    }

    public int getRPB() {
        return RPB;
    }

    public int getRPF() {
        return RPF;
    }

    public int getAPB() {
        return APB;
    }

    public int getAPF() {
        return APF;
    }

    public long getRawPosF() {
        return rawPosF;
    }

    public long getRawPosB() {
        return rawPosB;
    }

    public String getChrF() {
        return chrF;
    }

    public String getChrB() {
        return chrB;
    }

    public long getNumChrF() {
        return numChrF;
    }

    public long getNumChrB() {
        return numChrB;
    }

    public byte getStrandF() {
        return strandF;
    }

    public byte getStrandB() {
        return strandB;
    }

    public void setSvTypeCode(byte svTypeCode) {
        this.svTypeCode = svTypeCode;
    }

    public byte getSvTypeCode() {
        return svTypeCode;
    }

    public boolean isPpFlag() {
        return ppFlag;
    }

    public boolean isMmFlag() {
        return mmFlag;
    }

    public boolean isPmFlag() {
        return pmFlag;
    }

    public boolean isMpFlag() {
        return mpFlag;
    }

    public long getFrontCode() {
        return frontCode;
    }

    public long getBackCode() {
        return backCode;
    }
     
    public ArrayList<VariationV2> getVarList() {
        return varList;
    }   

    public Annotation getAnnoF() {
        return annoF;
    }

    public void setAnnoF(Annotation annoF) {
        this.annoF = annoF;
    }

    public Annotation getAnnoB() {
        return annoB;
    }

    public void setAnnoB(Annotation annoB) {
        this.annoB = annoB;
    }

    public boolean isSameChrFlag() {
        return sameChrFlag;
    }

    public boolean isDiffChrFlag() {
        return diffChrFlag;
    }

    public ArrayList<Double> getSameChrEnclideanDistance() {
        return sameChrEnclideanDistance;
    }

    public ArrayList<Double> getDiffChrEuclideanDistance() {
        return diffChrEuclideanDistance;
    }
    
    public void addEuclideanDistance(double euclideanDistanceValue){
        if(this.sameChrFlag==true){
            this.sameChrEnclideanDistance.add(euclideanDistanceValue);
        }else if(this.diffChrFlag == true){
            this.diffChrEuclideanDistance.add(euclideanDistanceValue);
        }
    }
    
    public void euclideanDistanceComplete(double euclideanDistanceValue){
        if(this.sameChrFlag==true){
            this.sameChrEnclideanDistanceCompleteFlag = true;
        }else if(this.diffChrFlag == true){
            this.diffChrEuclideanDistanceCompleteFlag = true;
        }
    }

    public long getEuCalFrontCode() {
        return euCalFrontCode;
    }

    public long getEuCalBackCode() {
        return euCalBackCode;
    }

    public boolean isTandemFlag() {
        return tandemFlag;
    }

    public void setTandemFlag(boolean tandemFlag) {
        this.tandemFlag = tandemFlag;
    }

    public boolean isDeleteFlag() {
        return deleteFlag;
    }

    public void setDeleteFlag(boolean deleteFlag) {
        this.deleteFlag = deleteFlag;
    }

    public boolean isIntraTransFlag() {
        return intraTransFlag;
    }

    public void setIntraTransFlag(boolean intraTransFlag) {
        this.intraTransFlag = intraTransFlag;
    }

    public boolean isInterTransFlag() {
        return interTransFlag;
    }

    public void setInterTransFlag(boolean interTransFlag) {
        this.interTransFlag = interTransFlag;
    }

    public boolean isChimericFlag() {
        return chimericFlag;
    }

    public void setChimericFlag(boolean chimericFlag) {
        this.chimericFlag = chimericFlag;
    }

    public boolean isAlienFlag() {
        return alienFlag;
    }

    public void setAlienFlag(boolean alienFlag) {
        this.alienFlag = alienFlag;
    }
    
    @Override
    public String toString(){
        return rawPosF+":"+strandF+"\t"+rawPosB+":"+strandB+"\t"+chrF+":"+RPF+"\t"+chrB+":"+RPB+"\t"+getNumCoverage()+"\t"+this.svTypeCode+":"+this.svType+"\t"+this.numStrandPattern;
    }
    
    public String shortSummary(){
        return chrF+":"+RPF+":"+strandF+"\t"+chrB+":"+RPB+":"+strandB+"\t"+getNumCoverage();
    }
  
}
