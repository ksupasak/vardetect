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
public class LASTResult {
    
    int score;
    String EG2;
    String E;

    String refChr;
    long refPos;
    int refNumBase;
    String refStrand;
    int refSize;
    String refSequenceMatch;
    long refSequenceCode;

    String sampleName;
    long samplePos;
    int sampleNumBase;
    String sampleStrand;
    int sampleSize;
    String sampleSequenceMatch;
    long sampleSequenceCode;
    
    public  LASTResult(){
        
    }
    
    public void addSampleData(String sampleName,long samplePos,int sampleNumBase,String sampleStrand,int sampleSize,String sampleSequenceMatch,long sampleSequenceCode){
        this.sampleName=sampleName;
        this.samplePos=samplePos;
        this.sampleNumBase=sampleNumBase;
        this.sampleStrand=sampleStrand;
        this.sampleSize=sampleSize;
        this.sampleSequenceMatch=sampleSequenceMatch;
        this.sampleSequenceCode=sampleSequenceCode;
    }
    
    public void addRefData(String refChr,long refPos,int refNumBase,String refStrand,int refSize,String refSequenceMatch,long refSequenceCode){
        this.refChr=refChr;
        this.refPos=refPos;
        this.refNumBase=refNumBase;
        this.refStrand=refStrand;
        this.refSize=refSize;
        this.refSequenceMatch=refSequenceMatch;
        this.refSequenceCode=refSequenceCode;
    }
    
    public void addScoreData(int score,String EG2,String E){
        this.score=score;
        this.EG2=EG2;
        this.E=E;
    }
    
    public void addKmer(int kmer){
        /**
        *   kmer is number of base would be use to decode per time
        */
    }
    
    @Override
    public String toString(){
        return "a "+"score="+this.score+" EG2="+this.EG2+" E="+this.E+"\n"
                + "s "+this.refChr+" "+this.refPos+" "+this.refNumBase+" "+this.refStrand+" "+this.refSize+" "+this.refSequenceMatch+"\n"
                + "s "+this.sampleName+" "+this.samplePos+" "+this.sampleNumBase+" "+this.sampleStrand+" "+this.sampleSize+" "+this.sampleSequenceMatch+"\n";
    }

    public int getScore() {
        return score;
    }

    public String getEG2() {
        return EG2;
    }

    public String getE() {
        return E;
    }

    public String getRefChr() {
        return refChr;
    }

    public long getRefPos() {
        return refPos;
    }

    public int getRefNumBase() {
        return refNumBase;
    }

    public String getRefStrand() {
        return refStrand;
    }

    public int getRefSize() {
        return refSize;
    }

    public String getRefSequenceMatch() {
        return refSequenceMatch;
    }

    public long getRefSequenceCode() {
        return refSequenceCode;
    }

    public String getSampleName() {
        return sampleName;
    }

    public long getSamplePos() {
        return samplePos;
    }

    public int getSampleNumBase() {
        return sampleNumBase;
    }

    public String getSampleStrand() {
        return sampleStrand;
    }

    public int getSampleSize() {
        return sampleSize;
    }

    public String getSampleSequenceMatch() {
        return sampleSequenceMatch;
    }

    public long getSampleSequenceCode() {
        return sampleSequenceCode;
    }
}
