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
public class SVGroup implements Comparable<SVGroup> {
    /**
     * Object that use to store group of object variationV2 which has the same structure variant pattern
     */
    private ArrayList<VariationV2> varList;
    private String svType;
    private int RPB;
    private int RPF;
    private int APB;
    private int APF;
    private long rawPosF;
    private long rawPosB;
    private String chrF;
    private String chrB;
    private byte strandF;
    private byte strandB;
    private int numCoverage;
    private byte svTypeCode;    // 0 is tanden, 1 is indel, 2 is intraTrans, 3 is interTrans, 4 is unclassify
    
    public SVGroup(){
        varList = new ArrayList();
    }
    
    public void addVariationV2(VariationV2 inVar){
        if(!this.varList.contains(inVar)){
            this.varList.add(inVar);
        } 
    }
    
    public String getSVType(){
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
        /**
         * Classify SV type
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

    public int getNumCoverage() {
        numCoverage = this.varList.size();
        return numCoverage;
    }

    @Override
    public int compareTo(SVGroup compareSVGroup) {
        int compareCov = ((SVGroup)compareSVGroup).getNumCoverage();
        /* For Ascending order*/
//        return this.studentage-compareage;
        /* For Descending order do like this */
        return compareCov-getNumCoverage();
    }

    public ArrayList<VariationV2> getVarList() {
        return varList;
    }   
    
    @Override
    public String toString(){
        return rawPosF+":"+strandF+"\t"+rawPosB+":"+strandB+"\t"+chrF+":"+RPF+"\t"+chrB+":"+RPB+"\t"+getNumCoverage()+"\t"+this.svTypeCode+":"+this.svType;
    }
  
}
