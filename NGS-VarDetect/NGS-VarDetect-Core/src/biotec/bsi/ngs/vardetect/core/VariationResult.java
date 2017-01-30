/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


/**
 *
 * @author worawich
 */
public class VariationResult {
    private Map<Integer,ArrayList<String[]>> variationResult;
    int numChrF; 
    long iniPosF; 
    long lastPosF;
    int iniIndexF;
    int greenF;
    int yellowF;
    int orangeF; 
    int redF; 
    int snpFlagF;
    int iniBackFlagF;
    String readNameF;
    String strandF;

    int numChrB; 
    long iniPosB; 
    long lastPosB;
    int iniIndexB;
    int greenB;
    int yellowB;
    int orangeB; 
    int redB; 
    int snpFlagB;
    int iniBackFlagB;
    String readNameB;
    String strandB;
    
    int merLength;
    int readLength;
    
    private ArrayList<Variation> listFusion;
    private ArrayList<Variation> listSNP;
    private ArrayList<Variation> listIndel;
    private ArrayList<Variation> listOthers;
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapFusion;
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapIndel;
    private Map<Long,ArrayList<Variation>> coverageMapSNP;
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapOther;
    
    private long mask28bit = 268435455;
    
    public VariationResult(){
        this.variationResult = new TreeMap();
        this.listFusion = new ArrayList();
        this.listSNP = new ArrayList();
        this.listIndel = new ArrayList();
        this.listOthers = new ArrayList();
        this.coverageMapFusion = new TreeMap();
        this.coverageMapIndel = new TreeMap();
        this.coverageMapOther = new TreeMap();
        this.coverageMapSNP = new TreeMap();
    }
    
    public void addMerLength(int merLen){
        this.merLength = merLen;
    }
    
    public void addReadLength(int readLen){
        this.readLength = readLen;
    }
    
    public void addVariationMap(Map<Integer,ArrayList<String[]>> variation){
//       this.variationResult.putAll(variation);
       
        Map<Integer,ArrayList<String[]>> firstMap = variation;
        Map<Integer,ArrayList<String[]>> secondMap = this.variationResult;
        
        Set<Map.Entry<Integer,ArrayList<String[]>>> entries = firstMap.entrySet();
        for ( Map.Entry<Integer,ArrayList<String[]>> entry : entries ) {
          ArrayList<String[]> secondMapValue = secondMap.get( entry.getKey() );
          if ( secondMapValue == null ) {
            secondMap.put( entry.getKey(), entry.getValue() );
          }
          else {
            secondMapValue.addAll( entry.getValue() );
          }
        }
        
        this.variationResult.putAll(secondMap);
  
    }
    
    public void createVariantReport(){
        this.variationResult.keySet();
        for(Map.Entry<Integer, ArrayList<String[]>> entry : this.variationResult.entrySet()) {
            /**
             * Loop over each type of variable00000
             */
            Integer key = entry.getKey();
            ArrayList<String[]> value = entry.getValue();
            for(int i=0;i<value.size();i++){
                String[] pairedPeak = value.get(i);
                String frontPeak = pairedPeak[0];
                String backPeak = pairedPeak[1];

                if(frontPeak != null){
                    String[] data = frontPeak.split(",");
                    numChrF = Integer.parseInt(data[0]);
                    iniPosF = Long.parseLong(data[1]);
                    lastPosF = Long.parseLong(data[2]);
                    iniIndexF = Integer.parseInt(data[8]);
                    greenF = Integer.parseInt(data[3]);
                    yellowF = Integer.parseInt(data[4]);
                    orangeF = Integer.parseInt(data[5]);
                    redF = Integer.parseInt(data[6]);
                    snpFlagF = Integer.parseInt(data[10]);
                    readNameF = data[9];
                    strandF = data[7];
                    iniBackFlagF = Integer.parseInt(data[11]);
                }
                if(backPeak != null){
                    String[] data = backPeak.split(",");
                    numChrB = Integer.parseInt(data[0]);
                    iniPosB = Long.parseLong(data[1]);
                    lastPosB = Long.parseLong(data[2]);
                    iniIndexB = Integer.parseInt(data[8]);
                    greenB = Integer.parseInt(data[3]);
                    yellowB = Integer.parseInt(data[4]);
                    orangeB = Integer.parseInt(data[5]);
                    redB = Integer.parseInt(data[6]);
                    snpFlagB = Integer.parseInt(data[10]);
                    readNameB = data[9];
                    strandB = data[7];
                    iniBackFlagB = Integer.parseInt(data[11]);
                }
                
                if(key == 0){
                    /**
                     * Variation type : SNP or other missing btw peak
                     */
                    Variation newVar = new Variation(this.merLength,this.readLength);
                    newVar.addType('S');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF);
                    this.listSNP.add(newVar);
                }else if(key == 1){
                    /**
                     * Variation type : Fusion
                     */
                    Variation newVar = new Variation(this.merLength,this.readLength);
                    newVar.addType('F');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB);
                    this.listFusion.add(newVar);
                }else if(key == 2){
                    /**
                     * Variation type : Indel both large and small
                     */
                    Variation newVar = new Variation(this.merLength,this.readLength);
                    newVar.addType('I');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB);
                    this.listIndel.add(newVar);
                }else if(key == 3){
                    /**
                     * Variation type : others with No green in both side (wasted)
                     */
                    Variation newVar = new Variation(this.merLength,this.readLength);
                    newVar.addType('O');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB);
                    this.listOthers.add(newVar);
                }
   
            }
        }
        
    }
    
    public void writeVarianReportToFile(String path , String nameFile) throws IOException{
         /**
         * Suitable for version 3 data structure (data structure that has iniIdx in its)
         * write result to file format for variant report
         */
        
        String filename = path+nameFile+".txt";
        PrintStream ps;
        FileWriter writer;        
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(filename,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        
        
        /**
         * Begin extract data to write
         */

        this.variationResult.keySet();
        for(Map.Entry<Integer, ArrayList<String[]>> entry : this.variationResult.entrySet()) {
            Integer key = entry.getKey();
            ArrayList<String[]> value = entry.getValue();
            
            if(key == 0){
                /**
                 * Variation type : SNP or other missing btw peak
                 */
                writer.write("Variation type : SNP or other missing btw peak");
                writer.write("\n"); 
            }else if(key == 1){
                /**
                 * Variation type : Fusion
                 */
                writer.write("Variation type : Fusion");
                writer.write("\n");
            }else if(key == 2){
                /**
                 * Variation type : Indel both large and small
                 */
                writer.write("Variation type : Indel both large and small");
                writer.write("\n");
            }else if(key == 3){
                /**
                 * Variation type : Indel both large and small
                 */
                writer.write("Variation type : Other case(wasted)");
                writer.write("\n");
            }
            
            for(int i=0;i<value.size();i++){
                String[] pairedPeak = value.get(i);
                String frontPeak = pairedPeak[0];
                String backPeak = pairedPeak[1];
                
                if(frontPeak != null){
                    String[] data = frontPeak.split(",");
                    numChrF = Integer.parseInt(data[0]);
                    iniPosF = Long.parseLong(data[1]);
                    lastPosF = Long.parseLong(data[2]);
                    iniIndexF = Integer.parseInt(data[8]);
                    greenF = Integer.parseInt(data[3]);
                    yellowF = Integer.parseInt(data[4]);
                    orangeF = Integer.parseInt(data[5]);
                    redF = Integer.parseInt(data[6]);
                    snpFlagF = Integer.parseInt(data[10]);
                    readNameF = data[9];
                    strandF = data[7];  
                }
                if(backPeak != null){
                    String[] data = backPeak.split(",");
                    numChrB = Integer.parseInt(data[0]);
                    iniPosB = Long.parseLong(data[1]);
                    lastPosB = Long.parseLong(data[2]);
                    iniIndexB = Integer.parseInt(data[8]);
                    greenB = Integer.parseInt(data[3]);
                    yellowB = Integer.parseInt(data[4]);
                    orangeB = Integer.parseInt(data[5]);
                    redB = Integer.parseInt(data[6]);
                    snpFlagB = Integer.parseInt(data[10]);
                    readNameB = data[9];
                    strandB = data[7]; 
                }
                writer.write("Variation"+i+" : ");
                writer.write(frontPeak + " || " + backPeak);
                writer.write("\n");
            }  
        }
 
        writer.flush();
        writer.close();
    }
    
    public void analyzeCoverageFusion(){
 
        for(int i =0;i<this.listFusion.size();i++){                         // loop over list of variation (Fusion type)
            Variation dummyVar = this.listFusion.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointFront();
            
            long bpFCode = (dummyVar.numChrF<<28)+bpF;
            long bpBCode = (dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapFusion.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapFusionII = this.coverageMapFusion.get(bpFCode);
                if(coverageMapFusionII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapFusionII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                }
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapFusionII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                coverageMapFusionII.put(bpBCode, coverageList);
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
            } 
        }
    }
    
    public void analyzeCoverageSNP(){
 
        for(int i =0;i<this.listSNP.size();i++){                         // loop over list of variation (SNP type)
            Variation dummyVar = this.listSNP.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointFront();
            
            long bpFCode = (dummyVar.numChrF<<28)+bpF;
//            long bpBCode = (dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapSNP.containsKey(bpFCode)){                      // check similarity of front breakpoint
                    
                /**
                 * put variation data to existing ArrayList<Variation>
                 * and put back to coverageMapFusionII
                 */

                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpFCode);
                coverageList.add(dummyVar);
                this.coverageMapSNP.put(bpFCode, coverageList);
 
            }else{
               
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                this.coverageMapSNP.put(bpFCode, coverageList);
                
            } 
        }
    }
    
    public void analyzeCoverageIndel(){
 
        for(int i =0;i<this.listIndel.size();i++){                         // loop over list of variation (Indel type)
            Variation dummyVar = this.listIndel.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointFront();
            
            long bpFCode = (dummyVar.numChrF<<28)+bpF;
            long bpBCode = (dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapIndel.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                if(coverageMapII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                }
                this.coverageMapIndel.put(bpFCode, coverageMapII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                coverageMapII.put(bpBCode, coverageList);
                this.coverageMapIndel.put(bpFCode, coverageMapII);
            } 
        }
    }
    
    public void analyzeCoverageOthers(){
 
        for(int i =0;i<this.listOthers.size();i++){                         // loop over list of variation (Others type)
            Variation dummyVar = this.listOthers.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointFront();
            
            long bpFCode = (dummyVar.numChrF<<28)+bpF;
            long bpBCode = (dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapOther.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapOther.get(bpFCode);
                if(coverageMapII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                }
                this.coverageMapOther.put(bpFCode, coverageMapII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                coverageMapII.put(bpBCode, coverageList);
                this.coverageMapOther.put(bpFCode, coverageMapII);
            } 
        }
    }
    
    public void writeVarianCoverageReportToFile(String path , String nameFile , char varType) throws IOException{
         /**
         * Suitable for version 3 data structure (data structure that has iniIdx in its)
         * write result to file format for variant report
         */
        
        String filename = path+nameFile+".txt";
        PrintStream ps;
        FileWriter writer;        
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(filename,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        
        int count = 1;
        if(varType == 'F'){
            writer.write("Variation type : Fusion\n");
            Set set = this.coverageMapFusion.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)bpFCode>>28;
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)bpBCode>>28;
                    
                    /**
                     * Write Report Part
                     */
                    writer.write("Group "+count);
                    writer.write("\tFront Break point : " + chrF +","+bpF);
                    writer.write("\tBack Break point : " + chrB +","+bpB);
 
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");
                    for(int i=0;i<coverageList.size();i++){
                        Variation var = coverageList.get(i);
                        
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                }
                count++;
            }
        }
        if(varType == 'I'){
            
        }
        if(varType == 'S'){
            
        }
        if(varType == 'O'){
            
        }

        writer.flush();
        writer.close();
    }
    
    
    
 
            
    
}
