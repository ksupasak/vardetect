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
public class ClusterGroup {
    
    private ArrayList<ShortgunSequence> clusterRead;
    private ArrayList<String> readNameList;                             // use for local alignment and other stuff
    private ArrayList<Byte> listChr;
    private ArrayList<Long> listIniPos;
    private ArrayList<Long> listLastPos;
    private ArrayList<Byte> listNumG;
    private ArrayList<Byte> listNumY;
    private ArrayList<Byte> listNumO;
    private ArrayList<Byte> listNumR;       
    private ArrayList<String> listStrand;       
    private ArrayList<Byte> listIniIndex;
    private int numMember;
           
            
    public ClusterGroup(){
        this.clusterRead = new ArrayList();
        this.readNameList = new ArrayList();
        this.listChr = new ArrayList();
        this.listIniIndex = new ArrayList();
        this.listIniPos = new ArrayList();
        this.listLastPos = new ArrayList();
        this.listNumG = new ArrayList();
        this.listNumY = new ArrayList();
        this.listNumO = new ArrayList();
        this.listNumR = new ArrayList();
        this.listStrand = new ArrayList();
    }
    
    public void addShortgunRead(ShortgunSequence readSS){
        this.clusterRead.add(readSS);
    }
    
    public ArrayList<ShortgunSequence> getShortgunRead(){
        return this.clusterRead;
    }
    
    public void addReadName(String inName){
        this.readNameList.add(inName);
    }
    
    public void addChromosomeNumber(Byte input){
        this.listChr.add(input);
    }
    
    public void addIniIndex(Byte input){
        this.listIniIndex.add(input);
    }
    
    public void addIniPos(long input){
        this.listIniPos.add(input);
    }
    
    public void addLastPos(long input){
        this.listLastPos.add(input);
    }
    
    public void addNumGreen(Byte input){
        this.listNumG.add(input);
    }
    
    public void addNumYellow(Byte input){
        this.listNumY.add(input);
    }
    
    public void addNumOrange(Byte input){
        this.listNumO.add(input);
    }
    
    public void addNumRed(Byte input){
        this.listNumR.add(input);
    }
    
    public void addStrand(String input){
        this.listStrand.add(input);
    }
    
    public int getNumMember(){
        this.numMember = this.readNameList.size();
        return this.numMember;
    }
    
    public ArrayList<Byte> getListChromosome(){
        return this.listChr;
    }
    
    public ArrayList<Long> getListIniPos(){
        return this.listIniPos;
    }
    
    public ArrayList<Long> getListLastPos(){
        return this.listLastPos;
    }
    
    public ArrayList<Byte> getListNumGreen(){
        return this.listNumG;
    }
    
    public ArrayList<Byte> getListNumYellow(){
        return this.listNumY;
    }
    
    public ArrayList<Byte> getListNumOrange(){
        return this.listNumO;
    }
    
    public ArrayList<Byte> getListNumRed(){
        return this.listNumR;
    }
    
    public ArrayList<String> getListStrand(){
        return this.listStrand;
    }
    
    public ArrayList<Byte> getListIniIndex(){
        return this.listIniIndex;
    }
    
    public ArrayList<String> getListReadname(){
        return this.readNameList;
    }
//    public void addGroupCharacteristic(long inChr, ){
//        
//    }
}
