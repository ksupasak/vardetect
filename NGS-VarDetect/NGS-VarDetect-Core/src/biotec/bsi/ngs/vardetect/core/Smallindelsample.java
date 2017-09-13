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
public class Smallindelsample {                                                       // Act as a storage for concatenate cut for easy get data
    private CharSequence sequence,cutA,cutB;
    private String type,chrA,chrB;
    private int iniA,iniB,indelSize,breakPointF,breakPointB; // initial position on each chrmosome
    private char indelType;
    
    public Smallindelsample(){
        
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
    public void addBasicInfo(String chrA, String chrB, int iniA, int iniB, char indelType, int indelSize){
        this.chrA = chrA;
        this.chrB = chrB;
        this.iniA = iniA;
        this.iniB = iniB;
        this.indelType = indelType;
        this.indelSize= indelSize;
    }
    public void addCutInfo(CharSequence cutA){
        this.cutA = cutA;     
    }
    public void addIndelInfo(CharSequence cutB){
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
    public CharSequence getIndelSequence(){
        return this.cutB;
    }
    public char getIndelType(){
        return this.indelType;
    }
    public int getIndelSize(){
        return this.indelSize;
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

