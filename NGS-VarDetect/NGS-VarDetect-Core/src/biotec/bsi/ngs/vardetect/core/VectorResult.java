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
public class VectorResult {
    private ArrayList<Long> listCode;
    private ArrayList<Long> listChr;
    private ArrayList<Long> listPos;
    private long chrPosCode1;
    private long chrPosCode2;
    private double vectorMagnitude; 
    private int numCode;
    
    public VectorResult(){
       listCode = new ArrayList();
       chrPosCode1 = 0;
       chrPosCode2 = 0;
       numCode = 0 ;
    }
    
    public void addChrPosCode(long code1,long code2){
        chrPosCode1 = code1;
        chrPosCode2 = code2;
        
        
        
        vectorMagnitude = Math.sqrt(Math.pow(code1,2)+Math.pow(code2,2));
        numCode = 2;
    }
    
    public void addChrPosCode(long code1){
        chrPosCode1 = code1;
        vectorMagnitude = Math.sqrt(Math.pow(code1,2)+0);
        numCode = 1;
        
    }
    
    //public ArrayList<Long> get
    
}
