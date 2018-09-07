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
public class LargeInsertionSample {
    
    private CharSequence sequenceF,sequenceB,cutA_F,cutB_F,cutA_B,cutB_B;
    private String type,chrA_F,chrB_F,chrA_B,chrB_B,oldCutA_F,oldCutB_F,oldCutA_B,oldCutB_B,snpBase_F,snpBase_B;
    private int iniA_F,iniB_F,iniA_B,iniB_B,breakPointF_F,breakPointB_F,breakPointF_B,breakPointB_B,insertionSize=0,posSNP_F,posSNP_B; // initial position on each chrmosome
    
    public LargeInsertionSample(){
        
    }
    
    public int getInsertionSize() {
        return insertionSize;
    }

    public void setInsertionSize(int insertionSize) {
        this.insertionSize = insertionSize;
    }
    
    public void addSequenceF(CharSequence in){
        this.sequenceF = in;
    }
    public void addSequenceB(CharSequence in){
        this.sequenceB = in;
    }
    public void addType(int inType){
        if (inType == 0){
            this.type = "+++";
        }else if(inType == 1){
            this.type = "+-+";
        }else if(inType == 2){
            this.type = "-+-";
        }else if(inType == 3){
            this.type = "---";
        }
    }
    public void addBasicInfo(String chrA_F, String chrB_F, String chrA_B, String chrB_B, int iniA_F, int iniB_F, int iniA_B, int iniB_B ){
        this.chrA_F = chrA_F;
        this.chrB_F = chrB_F;
        this.iniA_F = iniA_F;
        this.iniB_F = iniB_F;
        
        this.chrA_B = chrA_B;
        this.chrB_B = chrB_B;
        this.iniA_B = iniA_B;
        this.iniB_B = iniB_B;
    }
    public void addCutInfo(CharSequence cutA_F,CharSequence cutB_F,CharSequence cutA_B,CharSequence cutB_B){
        this.cutA_F = cutA_F;
        this.cutB_F = cutB_F;
        
        this.cutA_B = cutA_B;
        this.cutB_B = cutB_B;
    }
    public CharSequence getSequenceF(){
        return this.sequenceF;
    }
    public CharSequence getSequenceB(){
        return this.sequenceB;
    }
    public String getType(){
        return this.type;
    }
    public String getChrA_F() {
        return chrA_F;
    }

    public String getChrB_F() {
        return chrB_F;
    }

    public String getChrB_B() {
        return chrB_B;
    }

    public int getIniA_F() {
        return iniA_F;
    }

    public int getIniB_F() {
        return iniB_F;
    }

    public int getIniA_B() {
        return iniA_B;
    }

    public int getIniB_B() {
        return iniB_B;
    }

    public CharSequence getCutA_F() {
        return cutA_F;
    }

    public CharSequence getCutB_F() {
        return cutB_F;
    }

    public CharSequence getCutA_B() {
        return cutA_B;
    }

    public CharSequence getCutB_B() {
        return cutB_B;
    }

    public void addBreakPoint(int breakPointF_F, int breakPointB_F, int breakPointF_B, int breakPointB_B) {
        this.breakPointF_F = breakPointF_F;
        this.breakPointB_F = breakPointB_F;
        this.breakPointF_B = breakPointF_B;
        this.breakPointB_B = breakPointB_B;
    }

    public void setSNPInfo(int posSNP_F, String snpBase_F, int posSNP_B, String snpBase_B) {
        this.posSNP_F = posSNP_F;
        this.snpBase_F = snpBase_F;
        
        this.posSNP_B = posSNP_B;
        this.snpBase_B = snpBase_B;
    }

    public String getChrA_B() {
        return chrA_B;
    }

    public String getOldCutA_F() {
        return oldCutA_F;
    }

    public String getOldCutB_F() {
        return oldCutB_F;
    }

    public String getOldCutA_B() {
        return oldCutA_B;
    }

    public String getOldCutB_B() {
        return oldCutB_B;
    }

    public String getSnpBase_F() {
        return snpBase_F;
    }

    public String getSnpBase_B() {
        return snpBase_B;
    }

    public int getBreakPointF_F() {
        return breakPointF_F;
    }

    public int getBreakPointB_F() {
        return breakPointB_F;
    }

    public int getBreakPointF_B() {
        return breakPointF_B;
    }

    public int getBreakPointB_B() {
        return breakPointB_B;
    }

    public int getPosSNP_F() {
        return posSNP_F;
    }

    public int getPosSNP_B() {
        return posSNP_B;
    }

    public void addOldCutSequence(String oldCutA_F, String oldCutB_F, String oldCutA_B, String oldCutB_B) {
        this.oldCutA_F = oldCutA_F;
        this.oldCutB_F = oldCutB_F;
        
        this.oldCutA_B = oldCutA_B;
        this.oldCutB_B = oldCutB_B;
    }
}
