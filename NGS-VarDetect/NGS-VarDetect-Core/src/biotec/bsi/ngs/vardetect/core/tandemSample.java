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
public class tandemSample {
    private CharSequence sequence,cutA,cutB;
    private String type,chrA,chrB,snpBase;
    private int iniA,iniB,lastA,lastB,tandemSize,breakPointF,breakPointB,snpPos; // initial position on each chrmosome
    
    public tandemSample(){
        
    }
    public void addSequence(CharSequence in){
        this.sequence = in;
    }
    public void addBasicInfo(String chrA, String chrB, int iniA, int iniB, int lastA, int lastB, int tandemSize){
        this.chrA = chrA;
        this.chrB = chrB;
        this.iniA = iniA;
        this.iniB = iniB;
        this.lastA = lastA;
        this.lastB = lastB;
        this.tandemSize = tandemSize;
    }
    public void addBreakPoint(int breakPointF, int breakPointB) {
        this.breakPointF = breakPointF;
        this.breakPointB = breakPointB;
    }
    public void addCutInfo(CharSequence cutA,CharSequence cutB){
        this.cutA = cutA;
        this.cutB = cutB;
    }
    public void addType(int inType){
        if (inType == 0){
            this.type = "++";
        }else if(inType == 1){
            this.type = "--";
        }
    }
    public CharSequence getSequence() {
        return sequence;
    }

    public CharSequence getCutA() {
        return cutA;
    }

    public CharSequence getCutB() {
        return cutB;
    }

    public String getType() {
        return type;
    }

    public String getChrA() {
        return chrA;
    }

    public String getChrB() {
        return chrB;
    }

    public int getIniA() {
        return iniA;
    }

    public int getIniB() {
        return iniB;
    }

    public int getTandemSize() {
        return tandemSize;
    }

    public int getBreakPointF() {
        return breakPointF;
    }

    public int getBreakPointB() {
        return breakPointB;
    }

    public String getSnpBase() {
        return snpBase;
    }

    public int getSnpPos() {
        return snpPos;
    }
}
