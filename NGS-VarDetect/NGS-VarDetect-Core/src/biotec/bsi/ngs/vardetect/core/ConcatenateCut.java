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
public class ConcatenateCut {                                                       // Act as a storage for concatenate cut for easy get data
    private CharSequence sequence,cutA,cutB;
    private String type,chrA,chrB;
    private int iniA,iniB,breakPointF,breakPointB; // initial position on each chrmosome
    
    public ConcatenateCut(){
        
    }
    
    public void addSequence(CharSequence in){
        this.sequence = in;
    }
    public void addType(int inType){
        if (inType == 0){
            this.type = "++";
        }else if(inType == 1){
            this.type = "+-";
        }else if(inType == 2){
            this.type = "-+";
        }else{
            this.type = "--";
        }
    }
    public void addBasicInfo(String chrA, String chrB, int iniA, int iniB){
        this.chrA = chrA;
        this.chrB = chrB;
        this.iniA = iniA;
        this.iniB = iniB;
    }
    public void addCutInfo(CharSequence cutA,CharSequence cutB){
        this.cutA = cutA;
        this.cutB = cutB;
    }
    public CharSequence getSequence(){
        return this.sequence;
    }
    public String getType(){
        return this.type;
    }
    public int getiniA(){
        return this.iniA;
    }
    public int getiniB(){
        return this.iniB;
    }
    public String getchrA(){
        return this.chrA;
    }
    public String getchrB(){
        return this.chrB;
    }
    public CharSequence getcutA(){
        return this.cutA;
    }
    public CharSequence getcutB(){
        return this.cutB;
    }

    public void addBreakPoint(int breakPointF, int breakPointB) {
        this.breakPointF = breakPointF;
        this.breakPointB = breakPointB;
    }

    public int getBreakPointF() {
        return breakPointF;
    }

    public int getBreakPointB() {
        return breakPointB;
    }
}
