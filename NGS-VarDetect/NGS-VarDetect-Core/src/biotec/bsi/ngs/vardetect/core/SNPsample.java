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
public class SNPsample {                                                       // Act as a storage for concatenate cut for easy get data
    private CharSequence sequence,cutA,cutB;
    private String type,chrA;
    private int iniA,snpPos; // initial position on each chrmosome
    private char base;
    
    public SNPsample(){
        
    }
    
    public void addSequence(CharSequence in){
        this.sequence = in;
    }
    public void addSNPcharecteristic(int inType, char base , int snpPos){
        if (inType == 0){
            this.type = "+";
        }else if(inType == 1){
            this.type = "-";
        }
        this.base = base;
        this.snpPos = snpPos;
        
    }
    public void addBasicInfo(String chrA, int iniA){
        this.chrA = chrA;        
        this.iniA = iniA;       
    }
    public void addCutInfo(CharSequence cutA){
        this.cutA = cutA;   
    }
    public CharSequence getSequence(){
        return this.sequence;
    }
    public String getType(){
        return this.type;
    }
    
    public char baseChange(){
        return this.base;
    }
    
    public int getPosBaseChange(){
        return this.snpPos;
    }
    
    public int getiniA(){
        return this.iniA;
    }
    
    public String getchrA(){
        return this.chrA;
    }
    
    public CharSequence getcutA(){
        return this.cutA;
    }
    
    public char getBaseChange(){
        return this.baseChange();
    }
}
