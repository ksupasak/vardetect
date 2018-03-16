/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.Vector;

/**
 *
 * @author worawich
 */
public class Annotation {
    
    private String geneName;
    private long chrName;
    private String source;
    private String feature;
    private long start;
    private long stop;
    private String score;
    private String strand;
    private String frame;
    private String attribute;
    
    public Annotation(){

    }
    
    //Vector<ChromosomeSequence> chrs ;
    public Annotation(long inChrName,String inSource,String inFeature,long inStart,long inStop,String inScore,String inStrand,String inFrame,String inAttribute){

//        this.geneName = inGeneName;
        this.chrName = inChrName;
        this.source = inSource;
        this.feature = inFeature;
        this.start = inStart;
        this.stop = inStop;
        this.score = inScore;
        this.strand = inStrand;
        this.frame = inFrame;
        this.attribute = inAttribute;
        
    }
    
    public long getChrName(){
        return chrName;
    }
    public String getGeneName(){
        return geneName;
    }
    public String getSource(){
        return source;
    }
    public String getFeature(){
        return feature;
    }
    public long getStartPos(){
        return start;
    }
    public long getStopPos(){
        return stop;
    }
    public String getScore(){
        return score;
    }
    public String getStrand(){
        return strand;
    }
    public String getFrame(){
        return frame;
    }
    public String getAttribute(){
        return attribute;
    }
   
    @Override
    public String toString(){
        
        String str = chrName+" "+source+" "+feature+" "+strand+" "+attribute;
        return str;
    }
    
}
