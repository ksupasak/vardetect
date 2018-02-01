/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.toMap;

/**
 *
 * @author worawich
 */
public class VariationResult {
    private ArrayList<VariationV2> varList;
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
    int readLengthF;
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
    int readLengthB;
    String readNameB;
    String strandB;
    
    int merLength;
    int readLength;
    
    byte percentMatch;

    private ArrayList<Variation> listFusion;
    private ArrayList<Variation> listSNP;
    private ArrayList<Variation> listIndel;
    private ArrayList<Variation> listOthers;
    private ArrayList<Variation> listOneTail;                                                // store variation object of one tail event 
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapFusion;
    private Map<Long,Map<Long,ArrayList<Integer>>> sumNumMatchCoverageMapFusion;          // store ArrayList of Integer [sumNumMatchFront,sumNumMatchBack]. Use in annotation process to calculate average back and front breakpoint position use as search value for annotated gene
    private Map<String,Map<String,ArrayList<Integer>>> sumNumMatchCoverageMap;          // store ArrayList of Integer [sumNumMatchFront,sumNumMatchBack]. Use in annotation process to calculate average back and front breakpoint position use as search value for annotated gene
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapIndel;
    private Map<String,Map<String,SVGroup>> coverageMap;                   // coverageMap for new version of algorithm (speed up Version) did not pass SV type classification it sum all possible SV
    private Map<Long,ArrayList<Variation>> coverageMapSNP;
    private Map<Long,Map<Long,ArrayList<Variation>>> coverageMapOther;
    private Map<Long,ArrayList<Integer>> oneTailFrontLinkIdx;                                      // map that store key = original frontbreakpoint of one tail   and value is index of variaiton object store in listOnetail
    private Map<Long,ArrayList<Integer>> oneTailBackLinkIdx;                                       // map that store key = original backbreakpoint of one tail   and value is index of variaiton object store in listOnetail
    private ArrayList<Integer> listIndexOfUsedOneTail;                                      // store index of one tail that has been add to some indel coverage (we will used this index to remove those onetail from listOneTail when we want to export only one tail that doesn't match to any index type) not effect coverage discover part. this variable is stand for unmatch one-tail export part
    private ArrayList<ArrayList<Long>> sortedCoverageArrayFusion;
    private ArrayList<ArrayList<Long>> sortedCoverageArrayIndel;
    private ArrayList<SVGroup> tandemList;
    private ArrayList<SVGroup> indelList;
    private ArrayList<SVGroup> intraTransList;
    private ArrayList<SVGroup> interTransList;
    
    private long mask28bit = 268435455;
    
    public byte getPercentMatch() {
        return percentMatch;
    }

    public void setPercentMatch(byte percentMatch) {
        this.percentMatch = percentMatch;
    }
    
    public VariationResult(){
        this.variationResult = new TreeMap();
        this.listFusion = new ArrayList();
        this.listSNP = new ArrayList();
        this.listIndel = new ArrayList();
        this.listOthers = new ArrayList();
        this.listOneTail = new ArrayList();                                     
        this.coverageMapFusion = new LinkedHashMap();
        this.coverageMapIndel = new LinkedHashMap();
        this.coverageMap = new LinkedHashMap();
        this.coverageMapOther = new TreeMap();
        this.coverageMapSNP = new TreeMap();
        this.sumNumMatchCoverageMapFusion = new TreeMap();
        this.oneTailFrontLinkIdx = new TreeMap();
        this.oneTailBackLinkIdx = new TreeMap();
        this.listIndexOfUsedOneTail = new ArrayList();
        this.sortedCoverageArrayFusion = new ArrayList();
        this.sortedCoverageArrayIndel = new ArrayList();
        this.tandemList = new ArrayList();
        this.indelList = new ArrayList();
        this.intraTransList = new ArrayList();
        this.interTransList = new ArrayList();
    }
    
    public void addMerLength(int merLen){
        this.merLength = merLen;
    }
    
    public void addReadLength(int readLen){
        this.readLength = readLen;
    }

    public void setVarList(ArrayList<VariationV2> varList) {
        this.varList = varList;
    }
    
    public void addVariationMap(Map<Integer,ArrayList<String[]>> variation){
//       this.variationResult.putAll(variation);
       /**
        * This function stand for store all variationMap that come from different set of read 
        * come in different time. So the secondMap act as dummy for the data center that store all variaitonMap.
        */
        Map<Integer,ArrayList<String[]>> firstMap = variation;
        Map<Integer,ArrayList<String[]>> secondMap = this.variationResult;          // make secondMap as this.variationResult
        
        Set<Map.Entry<Integer,ArrayList<String[]>>> entries = firstMap.entrySet();
        for ( Map.Entry<Integer,ArrayList<String[]>> entry : entries ) {
          ArrayList<String[]> secondMapValue = secondMap.get( entry.getKey() );     // pull ArrayList<String[]> of specific indel type out for add new set of variaition
          if ( secondMapValue == null ) {
              /**
               * In case that the indel type is first time exist
               */
            secondMap.put( entry.getKey(), entry.getValue() );                      
          }
          else {
              /**
               * In case that the indel type is already exist. Just add data.
               */
            secondMapValue.addAll( entry.getValue() );
          }
        }
        
        this.variationResult.putAll(secondMap);
  
    }
    
    public void createVariantReport(){
        int countIndexOneTail = 0;
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
                    readLengthF = Integer.parseInt(data[12]);
                    
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
                    readLengthB = Integer.parseInt(data[12]);
                }
                
                if(key == 0){
                    /**
                     * Variation type : SNP or other missing btw peak
                     */
                    Variation newVar = new Variation(this.merLength);
                    newVar.addType('S');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);
                    this.listSNP.add(newVar);
                }else if(key == 1){
                    /**
                     * Variation type : Fusion
                     */
                    
                    Variation newVar = new Variation(this.merLength);
                    newVar.addType('F');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB, readLengthB);
                    this.listFusion.add(newVar);
                }else if(key == 2){
                    /**
                     * Variation type : Indel both large and small
                     */
                    Variation newVar = new Variation(this.merLength);
                    newVar.addType('I');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB, readLengthB);
                    this.listIndel.add(newVar);
                }else if(key == 3){
                    /**
                     * Variation type : others with No green in both side (wasted)
                     */
                    Variation newVar = new Variation(this.merLength);
                    newVar.addType('O');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);
                    newVar.addBackPeak(numChrB, iniPosB, lastPosB, greenB, yellowB, orangeB, redB, strandB, iniIndexB, readNameB, snpFlagB, iniBackFlagB, readLengthB);
                    this.listOthers.add(newVar);
                }else if(key == 4){
                    /**
                     * Variation type : 1 tail event (Not sure it is front or back part but data will come in with front variable)
                     * Because we cannot exactly decide that this one tail event is front or back. So, we will add the same data in to both front and back peak
                     * Which will calculate both back breakpoint and front breakpoint
                     */
                    Variation newVar = new Variation(this.merLength);
                    newVar.addType('T');
                    newVar.addFrontPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);
                    newVar.addBackPeak(numChrF, iniPosF, lastPosF, greenF, yellowF, orangeF, redF, strandF, iniIndexF, readNameF, snpFlagF, iniBackFlagF, readLengthF);                    
                    this.listOneTail.add(newVar);
                    
                    /**
                     * create Map to store index of each variation object that link with breakpoint front and back code (chromosome number plus with breakpoint)
                     * key = break point code
                     * value = arrayList of index on this.listOneTail
                     */
                    long bpFCode = ((long)newVar.numChrF<<28)+newVar.getOriBreakPointF();
                    long bpBCode = ((long)newVar.numChrB<<28)+newVar.getOriBreakPointB();

                    if(this.oneTailFrontLinkIdx.containsKey(bpFCode)){
                        ArrayList<Integer> listIndex = this.oneTailFrontLinkIdx.get(bpFCode);
                        listIndex.add(countIndexOneTail);
                        this.oneTailFrontLinkIdx.put(bpFCode, listIndex);
                    }else{
                        ArrayList<Integer> listIndex = new ArrayList();
                        listIndex.add(countIndexOneTail);
                        this.oneTailFrontLinkIdx.put(bpFCode, listIndex);
                    }
                    
                    if(this.oneTailBackLinkIdx.containsKey(bpBCode)){
                        ArrayList<Integer> listIndex = this.oneTailBackLinkIdx.get(bpBCode);
                        listIndex.add(countIndexOneTail);
                        this.oneTailBackLinkIdx.put(bpBCode, listIndex);
                    }else{
                        ArrayList<Integer> listIndex = new ArrayList();
                        listIndex.add(countIndexOneTail);
                        this.oneTailBackLinkIdx.put(bpBCode, listIndex);
                    }

                    countIndexOneTail++;
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
    
    public void analyzeCoverage(){
        /**
         * This fuction use to group same SV event togather
         * Use with new algorithm (speed align Version by P'soup)
         * Add remove duplication in coverage group (code has been implement in SVGroup object)
         */
        
        for(int i=0;i<this.varList.size();i++){
            
            VariationV2 dummyVar = varList.get(i);
            int bpF = dummyVar.getBreakpointF();
            int bpB = dummyVar.getBreakpointB();
            byte strandF = dummyVar.getStrandF();
            byte strandB = dummyVar.getStrandB();
            String chrF = dummyVar.getChrF();
            String chrB = dummyVar.getChrB();
       
            if(strandF==1 && strandB==1){
                /**
                 * For strand -- we will reverse the front and back part into normal arrange like reference 
                 * Front become back and back become front 
                 * Switch just for clustering task (Does not switch when classify SV type)
                 */
                bpF = dummyVar.getBreakpointB();
                bpB = dummyVar.getBreakpointF();
                chrF = dummyVar.getChrB();
                chrB = dummyVar.getChrF();
            }
            
            String bpFCode = chrF+":"+bpF;
            String bpBCode = chrB+":"+bpB;         
            
            if(this.coverageMap.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<String,SVGroup> coverageMapII = this.coverageMap.get(bpFCode);
//                Map<String,ArrayList<Integer>> sumNumMatchCoverageMapII = this.sumNumMatchCoverageMap.get(bpFCode);
                if(coverageMapII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    SVGroup svGroup = coverageMapII.get(bpBCode);
                    svGroup.addVariationV2(dummyVar);
                    coverageMapII.put(bpBCode, svGroup);
                    
//                    ArrayList<Integer> sumNumMatchList = sumNumMatchCoverageMapII.get(bpBCode);           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
//                    int sumNumMatchOldF = sumNumMatchList.get(0);
//                    int sumNumMatchOldB = sumNumMatchList.get(1);
//                    int sumNumMatchNewF = sumNumMatchOldF+dummyVar.numMatchF;
//                    int sumNumMatchNewB = sumNumMatchOldB+dummyVar.numMatchB;
//                    sumNumMatchList.add(0, sumNumMatchNewF);
//                    sumNumMatchList.add(1, sumNumMatchNewB);
//                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }else{
                    SVGroup svGroup = new SVGroup();
                    svGroup.addVariationV2(dummyVar);
                    coverageMapII.put(bpBCode, svGroup);
                    
//                    ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
//                    sumNumMatchList.add(0, dummyVar.numMatchF);
//                    sumNumMatchList.add(1, dummyVar.numMatchB);
//                    sumNumMatchCoverageMap.put(bpBCode, sumNumMatchList);
                }
                this.coverageMap.put(bpFCode, coverageMapII);
//                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapII);
            }else{
                Map<String,SVGroup> coverageMapII = new TreeMap();
                SVGroup svGroup = new SVGroup();
                svGroup.addVariationV2(dummyVar);
                
//                Map<String,ArrayList<Integer>> sumNumMatchCoverageMapII = new TreeMap();
//                ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
//                sumNumMatchList.add(0, dummyVar.numMatchF);
//                sumNumMatchList.add(1, dummyVar.numMatchB);
//                sumNumMatchCoverageMapII.put(bpBCode, sumNumMatchList);
//                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapII);
                
//                /**
//                 * Add data
//                 */
                coverageMapII.put(bpBCode, svGroup);
                this.coverageMap.put(bpFCode, coverageMapII);                
            }            
        }
        System.out.println("");
    }
    
    public void analyzeCoverageFusion(){
        /**
         * analyze coverage of Fusion event 
         * 1. group same front and back break point 
         * 2. consider one tail pattern
         * 
         * "Caution : this function must be called after createVariantReport function" 
         */
        
        for(int i =0;i<this.listFusion.size();i++){                         // loop over list of variation (Fusion type)
            Variation dummyVar = this.listFusion.get(i);
            long bpF = dummyVar.getOriBreakPointF();
            long bpB = dummyVar.getOriBreakPointB();
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
            long bpBCode = ((long)dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapFusion.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapFusionII = this.coverageMapFusion.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                if(coverageMapFusionII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapFusionII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = sumNumMatchCoverageMapFusionII.get(bpBCode);           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    int sumNumMatchOldF = sumNumMatchList.get(0);
                    int sumNumMatchOldB = sumNumMatchList.get(1);
                    int sumNumMatchNewF = sumNumMatchOldF+dummyVar.numMatchF;
                    int sumNumMatchNewB = sumNumMatchOldB+dummyVar.numMatchB;
                    sumNumMatchList.add(0, sumNumMatchNewF);
                    sumNumMatchList.add(1, sumNumMatchNewB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    sumNumMatchList.add(0, dummyVar.numMatchF);
                    sumNumMatchList.add(1, dummyVar.numMatchB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapFusionII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                
                
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = new TreeMap();
                ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                sumNumMatchList.add(0, dummyVar.numMatchF);
                sumNumMatchList.add(1, dummyVar.numMatchB);
                sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
                
//                /**
//                 * One-tail Coverage Discover
//                 * Check for one tail that has possibility to have the same either front or back break point to the breakpoint of this indel group
//                 * and put that into the coverage of this indel group (common 2-tail variation)
//                 * 
//                 * For first approach, we will not store numMatch information of one tail in sumNumMatch.
//                 */
//                
//                if(this.oneTailFrontLinkIdx.containsKey(bpFCode)){
//                    /**
//                     * In case that oneTail event is on front side
//                     */
//                    ArrayList<Integer> listIndex = this.oneTailFrontLinkIdx.get(bpFCode);
//                    
//                    for(Integer oneTailIndex : listIndex){
//                        Variation oneTail = this.listOneTail.get(oneTailIndex);
//                        coverageList.add(oneTail);
//                        
//                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
//                            // add onetail index into arraylist of used index
//                            this.listIndexOfUsedOneTail.add(oneTailIndex);
//                        }
//                    }
//                }else if(this.oneTailBackLinkIdx.containsKey(bpBCode)){
//                    /**
//                     * In case that oneTail event is on back side
//                     */
//                    ArrayList<Integer> listIndex = this.oneTailBackLinkIdx.get(bpBCode);
//                    
//                    for(Integer oneTailIndex : listIndex){
//                        Variation oneTail = this.listOneTail.get(oneTailIndex);
//                        coverageList.add(oneTail);
//                        
//                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
//                            // add onetail index into arraylist of used index
//                            this.listIndexOfUsedOneTail.add(oneTailIndex);
//                        }
//                    }                 
//                }
                /******************************************************************************************/
                
                /**
                 * Add data
                 */
                coverageMapFusionII.put(bpBCode, coverageList);
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
            } 
        }
        
        /**
         * Discover one-tail coverage
         * Loop existing coverageMapindel that already group 2-tail events 
         * use bpFCode and bpBCode as signature variable for grouping
         */
        Set set = this.coverageMapFusion.keySet();
        Iterator iterKey = set.iterator();
        while(iterKey.hasNext()){
            long bpFCode = (long)iterKey.next();
            long bpF = bpFCode&this.mask28bit;
            int chrF = (int)(bpFCode>>28);

            Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
            Set setII = coverageMapII.keySet();
            Iterator iterKeyII = setII.iterator();
            while(iterKeyII.hasNext()){
                long bpBCode = (long)iterKeyII.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                
                /**
                 * One-tail Coverage Discover
                 * Check for one tail that has possibility to have the same either front or back break point to the breakpoint of this indel group
                 * and put that into the coverage list of this indel group (common 2-tail variation)
                 * 
                 * For first approach, we will not store numMatch information of one tail in sumNumMatch.
                 */
                
//                if(dummyVar.readNameF.equals("Read01SS01")){
//                    System.out.println();
//                }
                
                if(this.oneTailFrontLinkIdx.containsKey(bpFCode)){
                    /**
                     * In case that oneTail event is on front side
                     */
                    ArrayList<Integer> listIndex = this.oneTailFrontLinkIdx.get(bpFCode);
                    
                    for(Integer oneTailIndex : listIndex){
                        Variation oneTail = this.listOneTail.get(oneTailIndex);
                        
                        if(containsVariation(coverageList,oneTail) == false){
                            // check oneTail variation object are already contain in coverage list or not. if not we will add it to coverageList
                            coverageList.add(oneTail);
                        }
                        
                        
                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
                            // add onetail index into arraylist of used index
                            this.listIndexOfUsedOneTail.add(oneTailIndex);
                        }
                    }
                }
                
                if(this.oneTailBackLinkIdx.containsKey(bpBCode)){
                    /**
                     * In case that oneTail event is on back side
                     */
                    ArrayList<Integer> listIndex = this.oneTailBackLinkIdx.get(bpBCode);
                    
                    for(Integer oneTailIndex : listIndex){
                        Variation oneTail = this.listOneTail.get(oneTailIndex);
                        if(containsVariation(coverageList,oneTail) == false){
                            // check oneTail variation object are already contain in coverage list or not. if not we will add it to coverageList
                            coverageList.add(oneTail);
                        }
                        
                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
                            // add onetail index into arraylist of used index
                            this.listIndexOfUsedOneTail.add(oneTailIndex);
                        }
                    }                 
                }
                /******************************************************************************************/
                coverageMapII.put(bpBCode, coverageList);
            }

            this.coverageMapFusion.put(bpFCode, coverageMapII);           
        }
    }
    
    public void analyzeCoverageFusionV2(){
        /**
         * analyze coverage of Fusion event 
         * 1. group same front and back break point 
         * 2. sum group that has switch front and back break point together (EX groupA has fb 100 and bb 200  has coverage 10 and groupB has fb 200 and bb 100 has coverage 8 ==> sum civerage of this two group = 18 read)
         * 3. not consider one tail pattern
         * 
         * "Caution : this function must be called after createVariantReport function" 
         */
        
        for(int i =0;i<this.listFusion.size();i++){                         // loop over list of variation (Fusion type)
            Variation dummyVar = this.listFusion.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointBack();
            
            /**
             * Consider Signature true breakpoint front and back
             * The true signature should have front breakpoint value lower than back breakpoint value
             * In other word we consider the event in a view of normal strand like a view of reference 
             */
            if(bpB-bpF < 0){
                /**
                 * bpB has lower value than bpF
                 * So,reassign breakpoint value in to the right way that it should be
                 * (Just switch bpF to bpB and bpB to bpF)
                 */
                long dummybpF = bpF;
                long dummybpB = bpB;
                bpF = dummybpB;
                bpB = dummybpF;
            }
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
            long bpBCode = ((long)dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapFusion.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapFusionII = this.coverageMapFusion.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                if(coverageMapFusionII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapFusionII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = sumNumMatchCoverageMapFusionII.get(bpBCode);           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    int sumNumMatchOldF = sumNumMatchList.get(0);
                    int sumNumMatchOldB = sumNumMatchList.get(1);
                    int sumNumMatchNewF = sumNumMatchOldF+dummyVar.numMatchF;
                    int sumNumMatchNewB = sumNumMatchOldB+dummyVar.numMatchB;
                    sumNumMatchList.add(0, sumNumMatchNewF);
                    sumNumMatchList.add(1, sumNumMatchNewB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapFusionII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    sumNumMatchList.add(0, dummyVar.numMatchF);
                    sumNumMatchList.add(1, dummyVar.numMatchB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapFusionII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                
                
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = new TreeMap();
                ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                sumNumMatchList.add(0, dummyVar.numMatchF);
                sumNumMatchList.add(1, dummyVar.numMatchB);
                sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);                
                
                /**
                 * Add data
                 */
                coverageMapFusionII.put(bpBCode, coverageList);
                this.coverageMapFusion.put(bpFCode, coverageMapFusionII);
            } 
        }
    }
    
    public void analyzeCoverageSNP(){
 
        for(int i =0;i<this.listSNP.size();i++){                         // loop over list of variation (SNP type)
            Variation dummyVar = this.listSNP.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointBack();
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
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
        /**
         * analyze coverage of Fusion event 
         * 1. group same front and back break point 
         * 2. consider one tail pattern
         * 
         * "Caution : this function must be called after createVariantReport function" 
         */
        
        for(int i =0;i<this.listIndel.size();i++){                         // loop over list of variation (Indel type)
            
            Variation dummyVar = this.listIndel.get(i);
            /**
             * Case check for insertion => if number of insertion base has more than number of match base on back part
             * we will skip it
             */
            if(dummyVar.getIndelType().equals("insert")){
                if(dummyVar.indelBase >= dummyVar.numMatchB){
                    continue;
                }
            }
            /******************************/
//            long bpF = dummyVar.getBreakPointFront();
//            long bpB = dummyVar.getBreakPointBack();
            long bpF = dummyVar.getOriBreakPointF();
            long bpB = dummyVar.getOriBreakPointB();
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
            long bpBCode = ((long)dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapIndel.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                if(coverageMapII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = sumNumMatchCoverageMapFusionII.get(bpBCode);           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    int sumNumMatchOldF = sumNumMatchList.get(0);
                    int sumNumMatchOldB = sumNumMatchList.get(1);
                    int sumNumMatchNewF = sumNumMatchOldF+dummyVar.numMatchF;
                    int sumNumMatchNewB = sumNumMatchOldB+dummyVar.numMatchB;
                    sumNumMatchList.add(0, sumNumMatchNewF);
                    sumNumMatchList.add(1, sumNumMatchNewB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    sumNumMatchList.add(0, dummyVar.numMatchF);
                    sumNumMatchList.add(1, dummyVar.numMatchB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }
                this.coverageMapIndel.put(bpFCode, coverageMapII);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = new TreeMap();
                ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                sumNumMatchList.add(0, dummyVar.numMatchF);
                sumNumMatchList.add(1, dummyVar.numMatchB);
                sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
                
//                /**
//                 * One-tail Coverage Discover
//                 * Check for one tail that has possibility to have the same either front or back break point to the breakpoint of this indel group
//                 * and put that into the coverage of this indel group (common 2-tail variation)
//                 * 
//                 * For first approach, we will not store numMatch information of one tail in sumNumMatch.
//                 */
//                
////                if(dummyVar.readNameF.equals("Read01SS01")){
////                    System.out.println();
////                }
//                
//                if(this.oneTailFrontLinkIdx.containsKey(bpFCode)){
//                    /**
//                     * In case that oneTail event is on front side
//                     */
//                    ArrayList<Integer> listIndex = this.oneTailFrontLinkIdx.get(bpFCode);
//                    
//                    for(Integer oneTailIndex : listIndex){
//                        Variation oneTail = this.listOneTail.get(oneTailIndex);
//                        coverageList.add(oneTail);
//                        
//                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
//                            // add onetail index into arraylist of used index
//                            this.listIndexOfUsedOneTail.add(oneTailIndex);
//                        }
//                    }
//                }
//                
//                if(this.oneTailBackLinkIdx.containsKey(bpBCode)){
//                    /**
//                     * In case that oneTail event is on back side
//                     */
//                    ArrayList<Integer> listIndex = this.oneTailBackLinkIdx.get(bpBCode);
//                    
//                    for(Integer oneTailIndex : listIndex){
//                        Variation oneTail = this.listOneTail.get(oneTailIndex);
//                        coverageList.add(oneTail);
//                        
//                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
//                            // add onetail index into arraylist of used index
//                            this.listIndexOfUsedOneTail.add(oneTailIndex);
//                        }
//                    }                 
//                }
//                /******************************************************************************************/
//                
//                /**
//                 * Add data
//                 */
                coverageMapII.put(bpBCode, coverageList);
                this.coverageMapIndel.put(bpFCode, coverageMapII);
                
            } 
        }
        
        /**
         * Discover one-tail coverage
         * Loop existing coverageMapindel that already group 2-tail events 
         * use bpFCode and bpBCode as signature variable for grouping
         */
        Set set = this.coverageMapIndel.keySet();
        Iterator iterKey = set.iterator();
        while(iterKey.hasNext()){
            long bpFCode = (long)iterKey.next();
            long bpF = bpFCode&this.mask28bit;
            int chrF = (int)(bpFCode>>28);

            Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
            Set setII = coverageMapII.keySet();
            Iterator iterKeyII = setII.iterator();
            while(iterKeyII.hasNext()){
                long bpBCode = (long)iterKeyII.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                
                /**
                 * One-tail Coverage Discover
                 * Check for one tail that has possibility to have the same either front or back break point to the breakpoint of this indel group
                 * and put that into the coverage list of this indel group (common 2-tail variation)
                 * 
                 * For first approach, we will not store numMatch information of one tail in sumNumMatch.
                 */
                
//                if(dummyVar.readNameF.equals("Read01SS01")){
//                    System.out.println();
//                }
                
                if(this.oneTailFrontLinkIdx.containsKey(bpFCode)){
                    /**
                     * In case that oneTail event is on front side
                     */
                    ArrayList<Integer> listIndex = this.oneTailFrontLinkIdx.get(bpFCode);
                    
                    for(Integer oneTailIndex : listIndex){
                        Variation oneTail = this.listOneTail.get(oneTailIndex);
                        
                        if(containsVariation(coverageList,oneTail) == false){
                            // check oneTail variation object are already contain in coverage list or not. if not we will add it to coverageList
                            oneTail.setIndelType("front_oneTail");
                            coverageList.add(oneTail);
                        }
                        
                        
                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
                            // add onetail index into arraylist of used index
                            this.listIndexOfUsedOneTail.add(oneTailIndex);
                        }
                    }
                }
                
                if(this.oneTailBackLinkIdx.containsKey(bpBCode)){
                    /**
                     * In case that oneTail event is on back side
                     */
                    ArrayList<Integer> listIndex = this.oneTailBackLinkIdx.get(bpBCode);
                    
                    for(Integer oneTailIndex : listIndex){
                        Variation oneTail = this.listOneTail.get(oneTailIndex);
                        if(containsVariation(coverageList,oneTail) == false){
                            // check oneTail variation object are already contain in coverage list or not. if not we will add it to coverageList
                            oneTail.setIndelType("back_oneTail");
                            coverageList.add(oneTail);
                        }
                        
                        if(!this.listIndexOfUsedOneTail.contains(oneTailIndex)){        // if it not exist add new one
                            // add onetail index into arraylist of used index
                            this.listIndexOfUsedOneTail.add(oneTailIndex);
                        }
                    }                 
                }
                /******************************************************************************************/
                coverageMapII.put(bpBCode, coverageList);
            }

            this.coverageMapIndel.put(bpFCode, coverageMapII);           
        }
    }
    
    public void analyzeCoverageIndelV2(){
        
        /**
         * analyze coverage of Indel event 
         * 1. group same front and back break point 
         * 2. sum group that has switch front and back break point together (EX groupA has fb 100 and bb 200  has coverage 10 and groupB has fb 200 and bb 100 has coverage 8 ==> sum coverage of this two group = 18 read)
         * 3. not consider one tail pattern
         * 4. create parallel coverageMap that store the information of different strand pattern that present in the group of coverage 
         * 
         * "Caution : this function must be called after createVariantReport function"
         */
        
        for(int i =0;i<this.listIndel.size();i++){                         // loop over list of variation (Indel type)
            
            Variation dummyVar = this.listIndel.get(i);
            /**
             * Case check for insertion => if number of insertion base has more than number of match base on back part
             * we will skip it
             */
            if(dummyVar.getIndelType().equals("insert")){
                if(dummyVar.indelBase >= dummyVar.numMatchB){
                    continue;
                }
            }
            /******************************/
//            long bpF = dummyVar.getBreakPointFront();
//            long bpB = dummyVar.getBreakPointBack();
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointBack();
            String strandPattern = dummyVar.strandF + dummyVar.strandB;
            
            if(strandPattern.equals("++") || strandPattern.equals("--")){
                /**
                 * Consider only pattern that has strand ++ and -- (for +- and -+ it is different pattern use the original one [no switch])
                 * Consider Signature true breakpoint front and back
                 * The true signature should have front breakpoint value lower than back breakpoint value
                 * In other word we consider the event in a view of normal strand like a view of reference 
                 */
                if(bpB-bpF < 0){
                    /**
                     * bpB has lower value than bpF
                     * So,re-assign breakpoint value in to the right way that it should be (order from low to high)
                     * (Just switch bpF to bpB and bpB to bpF)
                     */
                    long dummybpF = bpF;
                    long dummybpB = bpB;
                    bpF = dummybpB;
                    bpB = dummybpF;
                }
            }
            
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
            long bpBCode = ((long)dummyVar.numChrB<<28)+bpB;
            
            if(this.coverageMapIndel.containsKey(bpFCode)){                      // check similarity of front breakpoint
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                if(coverageMapII.containsKey(bpBCode)){                       // check similarity of back breakpoit
                    /**
                     * put variation data to existing ArrayList<Variation>
                     * and put back to coverageMapFusionII
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = sumNumMatchCoverageMapFusionII.get(bpBCode);           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    int sumNumMatchOldF = sumNumMatchList.get(0);
                    int sumNumMatchOldB = sumNumMatchList.get(1);
                    int sumNumMatchNewF = sumNumMatchOldF+dummyVar.numMatchF;
                    int sumNumMatchNewB = sumNumMatchOldB+dummyVar.numMatchB;
                    sumNumMatchList.add(0, sumNumMatchNewF);
                    sumNumMatchList.add(1, sumNumMatchNewB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }else{
                    ArrayList<Variation> coverageList = new ArrayList();
                    coverageList.add(dummyVar);
                    coverageMapII.put(bpBCode, coverageList);
                    
                    ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                    sumNumMatchList.add(0, dummyVar.numMatchF);
                    sumNumMatchList.add(1, dummyVar.numMatchB);
                    sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                }
                this.coverageMapIndel.put(bpFCode, coverageMapII);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
            }else{
                Map<Long,ArrayList<Variation>> coverageMapII = new TreeMap();
                ArrayList<Variation> coverageList = new ArrayList();
                coverageList.add(dummyVar);
                
                Map<Long,ArrayList<Integer>> sumNumMatchCoverageMapFusionII = new TreeMap();
                ArrayList<Integer> sumNumMatchList = new ArrayList();           // ArrayList of integer is store 2 value sumNumMatchF (index 0)  and  sumNuMatchB (index 1)
                sumNumMatchList.add(0, dummyVar.numMatchF);
                sumNumMatchList.add(1, dummyVar.numMatchB);
                sumNumMatchCoverageMapFusionII.put(bpBCode, sumNumMatchList);
                this.sumNumMatchCoverageMapFusion.put(bpFCode, sumNumMatchCoverageMapFusionII);
                
//                /**
//                 * Add data
//                 */
                coverageMapII.put(bpBCode, coverageList);
                this.coverageMapIndel.put(bpFCode, coverageMapII);
                
            } 
        }
    }
    
    public void analyzeCoverageOthers(){
 
        for(int i =0;i<this.listOthers.size();i++){                         // loop over list of variation (Others type)
            Variation dummyVar = this.listOthers.get(i);
            long bpF = dummyVar.getBreakPointFront();
            long bpB = dummyVar.getBreakPointBack();
            
            long bpFCode = ((long)dummyVar.numChrF<<28)+bpF;
            long bpBCode = ((long)dummyVar.numChrB<<28)+bpB;
            
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
    
    public void writeVariantCoverageReportToFile(String nameFile , int coverageThreshold , char varType) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for variant report
        * 
        * Add option on forth argument : Boolean variable to turn on or off generate read name report of match 2-tail and 1-tail both un-match and match
        * 
        * "Caution : this function must be called after analyzeCoverage[Indel or Fusion] function"
        */
        
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_Fusion"+".txt";
                break;
            case 'I':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_Indel"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_SNP"+".txt";
                break;
            default:
                break;
        }

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
        
        
        
        if(varType == 'F'){
            int count = 1;
            writer.write("Variation type : Fusion\n");
            Set set = this.coverageMapFusion.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        writer.write("Group "+count);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB); 
    //                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write("\n");
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writer.write("\n");
                                }
                            }else{
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writer.write("\n");
                                
                            }
                        }
                        count++;
                    }
                }  
            }
        }
        if(varType == 'I'){
            int count = 1;
            writer.write("Variation type : Indel\n");
            Set set = this.coverageMapIndel.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        /**
                         * pick one variation to get indel information 
                         */
                        Variation tempVar = coverageList.get(0);
                        String indelType = tempVar.getIndelType();
                        long indelBase = tempVar.getIndelBase();
                        /***/
                        
                        writer.write("Group "+count);
                        writer.write("\tIndel Type : " + indelType);
                        writer.write("\tIndel Base : " + indelBase);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            
                            Variation var = coverageList.get(i);                            
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write("\n");
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writer.write("\n");
                                }
                                
                            }else{
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writer.write("\n");
                            }
                        }
                        count++;
                    }
                }
                
            }
            
        }
        if(varType == 'S'){
            int count = 1;
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +":"+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);

                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            
        }
        if(varType == 'O'){
            
        }

        writer.flush();
        writer.close();
    }
    
    public void writeVariantCoverageReportToFile(String nameFile , int coverageThreshold , char varType, boolean rnoReportFlag) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for variant report
        * 
        * Add option on forth argument : Boolean variable to turn on or off generate read name report of match 2-tail and 1-tail both un-match and match (rno : read name only report [report that contain only read name that has indel event])
        * 
        * "Caution : this function must be called after analyzeCoverage[Indel or Fusion] function"
        */
        ArrayList<String> usedReadName = new ArrayList();
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        String rnoFileName = "";                            // read name only report file name
        
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_Fusion"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_Fusion"+".txt";
                break;
            case 'I':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_Indel"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_Indel"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_CoverageReport_SNP"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_SNP"+".txt";
                break;
            default:
                break;
        }

        PrintStream ps;
        FileWriter writer;
        FileWriter rnoWriter;
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
        
        File rnoFile = new File(rnoFileName); //File object        
        if(rnoFile.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            rnoWriter = new FileWriter(rnoFileName,true);
        }else{
//            ps = new PrintStream(filename);
            rnoWriter = new FileWriter(rnoFileName);
        }
        
        
        
        if(varType == 'F'){
            int count = 1;
            writer.write("Variation type : Fusion\n");
            Set set = this.coverageMapFusion.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        writer.write("Group "+count);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB); 
    //                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write("\n");
                                    
                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                    
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writer.write("\n");
                                    
                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }
                            }else{
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writer.write("\n");
                                
                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                } 
                            }
                        }
                        count++;
                    }
                }  
            }
        }
        if(varType == 'I'){
            int count = 1;
            writer.write("Variation type : Indel\n");
            Set set = this.coverageMapIndel.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        /**
                         * pick one variation to get indel information 
                         */
                        Variation tempVar = coverageList.get(0);
                        String indelType = tempVar.getIndelType();
                        long indelBase = tempVar.getIndelBase();
                        /***/
                        
                        writer.write("Group "+count);
                        writer.write("\tIndel Type : " + indelType);
                        writer.write("\tIndel Base : " + indelBase);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            
                            Variation var = coverageList.get(i);

                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write("\n");
                                    
                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writer.write(" || ");
                                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writer.write("\n");
                                    
                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }
                                
                            }else{
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writer.write("\n");
                                
                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }
                            }
                        }
                        count++;
                    }
                }
                
            }
            
        }
        if(varType == 'S'){
            int count = 1;
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +":"+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);

                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            
        }
        if(varType == 'O'){
            
        }

        writer.flush();
        writer.close();
        
        rnoWriter.flush();
        rnoWriter.close();
    }
    
    public void writeVariantSortedCoverageReportToFile(String nameFile , int coverageThreshold , char varType, boolean rnoReportFlag) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for variant report in sorted order (high to low)
        * 
        * Add option on forth argument : Boolean variable to turn on or off generate read name report of match 2-tail and 1-tail both un-match and match (rno : read name only report [report that contain only read name that has indel event])
        * 
        * "Caution : this function must be called after sortCoverage[Indel or Fusion] function"
        */
        ArrayList<String> usedReadName = new ArrayList();
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        String rnoFileName = "";                            // read name only report file name
        String sdFile = "";
        String idFile = "";
        String ldFile = "";
        
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_SortedCoverageReport_Fusion"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_Fusion"+".txt";
                break;
            case 'I':
                sdFile = dummy[0]+"_match"+this.percentMatch+"_SortedCoverageReport_smallDelete"+".txt";
                idFile = dummy[0]+"_match"+this.percentMatch+"_SortedCoverageReport_smallInsert"+".txt";
                ldFile = dummy[0]+"_match"+this.percentMatch+"_SortedCoverageReport_largeDelete"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_Indel"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_SortedCoverageReport_SNP"+".txt";
                rnoFileName = dummy[0]+"_match"+this.percentMatch+"_readNameReport_SNP"+".txt";
                break;
            default:
                break;
        }

        PrintStream ps;
        FileWriter writer;
        FileWriter writerSD;    // small indel
        FileWriter writerID;    // small insertion
        FileWriter writerLD;    // large indel
        FileWriter rnoWriter;
        /**
         * Check File existing
         */

        File rnoFile = new File(rnoFileName); //File object        
        if(rnoFile.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            rnoWriter = new FileWriter(rnoFileName,true);
        }else{
//            ps = new PrintStream(filename);
            rnoWriter = new FileWriter(rnoFileName);
        }
        
        
        
        if(varType == 'F'){
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
            writer.write("Variation type : Fusion\n");
            for(int i=0;i<this.sortedCoverageArrayFusion.size();i++){
                if(this.sortedCoverageArrayFusion.size() == 0){
                    break;
                }
                ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayFusion.get(i);
                long bpFCode = dummyCoverageList.get(1);
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                long bpBCode = dummyCoverageList.get(2);
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                
                ArrayList<Variation> coverageList = getFusionCoverageList(bpFCode,bpBCode);
                
                if(coverageList.size()>coverageThreshold){
                    writer.write("Group "+count);
                    writer.write("\tFront Break point : " + chrF +":"+bpF);
                    writer.write("\tBack Break point : " + chrB +":"+bpB); 
//                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");
                    for(int j=0;j<coverageList.size();j++){
                        Variation var = coverageList.get(j);

                        if(var.variationType == 'T'){
                            /**
                            * Check one tail type
                            * to classify it side Front or Back
                            * and write report
                            */
                            if(var.getOriBreakPointF()==bpF){
                                /**
                                 * This one tail is sit on the front
                                 */
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                writer.write("\n");

                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }

                            }else if(var.getOriBreakPointB()==bpB){
                                /**
                                 * This one tail is sit on the back
                                 */
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                writer.write(" || ");
                                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writer.write("\n");

                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }
                            }
                        }else{
                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                            writer.write(" || ");
                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                            writer.write("\n");

                            if(!usedReadName.contains(var.readNameF)){
                                rnoWriter.write(var.readNameF);
                                rnoWriter.write("\n");
                                usedReadName.add(var.readNameF);
                            } 
                        }
                    }
                }
                count++;
            }
            writer.flush();
            writer.close();
        }
        if(varType == 'I'){
            writerSD = new FileWriter(sdFile);  //small delete File writer
            writerID = new FileWriter(idFile);  //small insert File writer
            writerLD = new FileWriter(ldFile);  //large delete File writer
            
            int countSD = 1;
            int countID = 1;
            int countLD = 1;

//            writer.write("Variation type : Indel\n");
            
            for(int i=0;i<this.sortedCoverageArrayIndel.size();i++){
                if(this.sortedCoverageArrayIndel.size() == 0){
                    break;
                }
                ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayIndel.get(i);
                long bpFCode = dummyCoverageList.get(1);
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                long bpBCode = dummyCoverageList.get(2);
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = getIndelCoverageList(bpFCode,bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                     * pick one variation to get indel information 
                     */
                    Variation tempVar = coverageList.get(0);
                    String indelType = tempVar.getIndelType();
                    long indelBase = tempVar.getIndelBase();
                    /***/
                    
                    if(indelType.equals("delete")){
                        writerSD.write("Group "+countSD);
                        writerSD.write("\tIndel Type : " + indelType);
                        writerSD.write("\tIndel Base : " + indelBase);
                        writerSD.write("\tFront Break point : " + chrF +":"+bpF);
                        writerSD.write("\tBack Break point : " + chrB +":"+bpB);
                        writerSD.write("\tCoverage : " + coverageList.size());
                        writerSD.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);
   
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writerSD.write(" || ");
                                    writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerSD.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerSD.write(" || ");
                                    writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writerSD.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }

                            }else{
                                writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writerSD.write(" || ");
                                writerSD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writerSD.write("\n");

                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }
                            }
                        }
                        countSD++;
                    }else if(indelType.equals("insert")){
                        writerID.write("Group "+countID);
                        writerID.write("\tIndel Type : " + indelType);
                        writerID.write("\tIndel Base : " + indelBase);
                        writerID.write("\tFront Break point : " + chrF +":"+bpF);
                        writerID.write("\tBack Break point : " + chrB +":"+bpB);
                        writerID.write("\tCoverage : " + coverageList.size());
                        writerID.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);
    
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writerID.write(" || ");
                                    writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerID.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerID.write(" || ");
                                    writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writerID.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }

                            }else{
                                writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writerID.write(" || ");
                                writerID.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writerID.write("\n");

                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }
                            }
                        }
                        countID++;
                    }else{
                        // large indel

                        writerLD.write("Group "+countLD);
                        writerLD.write("\tIndel Type : " + indelType);
                        writerLD.write("\tIndel Base : " + indelBase);
                        writerLD.write("\tFront Break point : " + chrF +":"+bpF);
                        writerLD.write("\tBack Break point : " + chrB +":"+bpB);
                        writerLD.write("\tCoverage : " + coverageList.size());
                        writerLD.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);
    
                            if(var.variationType == 'T'){
                                /**
                                * Check one tail type
                                * to classify it side Front or Back
                                * and write report
                                */
                                if(var.getOriBreakPointF()==bpF){
                                    /**
                                     * This one tail is sit on the front
                                     */
                                    writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                    writerLD.write(" || ");
                                    writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerLD.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }else if(var.getOriBreakPointB()==bpB){
                                    /**
                                     * This one tail is sit on the back
                                     */
                                    writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", 0,0,0,0,0,0,0,null,0,null,0,0,0));
                                    writerLD.write(" || ");
                                    writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                    writerLD.write("\n");

                                    if(!usedReadName.contains(var.readNameF)){
                                        rnoWriter.write(var.readNameF);
                                        rnoWriter.write("\n");
                                        usedReadName.add(var.readNameF);
                                    }
                                }

                            }else{
                                writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
                                writerLD.write(" || ");
                                writerLD.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
                                writerLD.write("\n");

                                if(!usedReadName.contains(var.readNameF)){
                                    rnoWriter.write(var.readNameF);
                                    rnoWriter.write("\n");
                                    usedReadName.add(var.readNameF);
                                }
                            }
                        }
                        countLD++;
                    }
                }  
            }
            writerSD.flush();
            writerSD.close();
            
            writerID.flush();
            writerID.close();
            
            writerLD.flush();
            writerLD.close();  
        }
        if(varType == 'S'){
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
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +":"+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);

                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            writer.flush();
            writer.close();
        }
        if(varType == 'O'){
            
        }


        
        rnoWriter.flush();
        rnoWriter.close();
    }
    
    public void writeVariantCoverageVirtualizeReportToFile(String nameFile , int coverageThreshold , char varType) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for variant report
        * create virtual sequence (give a point of view of support read)
        */
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageReport_Fusion"+".txt";
                break;
            case 'I':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageReport_Indel"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageReport_SNP"+".txt";
                break;
            default:
                break;
        }
        
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
        
        
        if(varType == 'F'){
            int count = 1;
            writer.write("Variation type : Fusion\n");
            Set set = this.coverageMapFusion.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        writer.write("Group "+count);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB); 
    //                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            writer.write(var.virtualSequence());
                            writer.write("\n");
                        }
                        count++;
                    }
                }  
            }
        }
        if(varType == 'I'){
            int count = 1;
            writer.write("Variation type : Indel\n");
            Set set = this.coverageMapIndel.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        /**
                         * pick one variation to get indel information 
                         */
                        Variation tempVar = coverageList.get(0);
                        String indelType = tempVar.getIndelType();
                        long indelBase = tempVar.getIndelBase();
                        /***/
                        
                        writer.write("Group "+count);
                        writer.write("\tIndel Type : " + indelType);
                        writer.write("\tIndel Base : " + indelBase);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            writer.write(var.virtualSequence());
                            writer.write("\n");
//                            Variation var = coverageList.get(i);
//                            
//                            if(var.variationType == 'T'){
//                                /**
//                                * Check one tail type
//                                * to classify it side Front or Back
//                                * and write report
//                                */
//                                if(var.getOriBreakPointF()==bpF){
//                                    /**
//                                     * This one tail is sit on the front
//                                     */
//                                    writer.write(var.oneTailVirtualSequenceFront());
//                                    writer.write("\n");
//                                }else if(var.getOriBreakPointB()==bpB){
//                                    /**
//                                     * This one tail is sit on the back
//                                     */
////                                    writer.write(var.oneTailVirtualSequenceBack());
//                                    writer.write("\n");
//                                }
//                            }else{
//                                writer.write(var.virtualSequence());
//                                writer.write("\n");
//                            }
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
//                            writer.write(" || ");
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
//                            writer.write("\n");
                        }
                        count++;
                    }
                }
                
            }
            
        }
        if(varType == 'S'){
            int count = 1;
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +","+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);
                        
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            
        }
        if(varType == 'O'){
            
        }

        writer.flush();
        writer.close();
    }
    
    public void writeVariantCoverageVirtualizeWithAnnotationReportToFile(String nameFile ,String gffFile, int coverageThreshold , char varType) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for variant report
        * create virtual sequence (give a point of view of support read)
        */
        
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV2(gffFile, "gene");
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageAnnotatedReport_Fusion"+".txt";
                break;
            case 'I':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageAnnotatedReport_Indel"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_VirtualCoverageAnnotatedReport_SNP"+".txt";
                break;
            default:
                break;
        }
        
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
        
        
        if(varType == 'F'){
            int count = 1;
            writer.write("Variation type : Fusion\n");
            Set set = this.coverageMapFusion.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchMapII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    ArrayList<Integer> sumNumMatchList = sumNumMatchMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        /**
                         * Calculate average numMatch for Front and Back part and Map to annotation binary tree for gene information
                         * create code chrPosStart and chrPosStop both front and back (use for binary search)
                         */
                        
                        int avgNumMatchF = sumNumMatchList.get(0)/coverageList.size();
                        int avgNumMatchB = sumNumMatchList.get(1)/coverageList.size();
                        
//                        long avgIniPosF = bpF-avgNumMatchF;
//                        long avgLastPosB = bpB+avgNumMatchB;
                        
                        long avgIniPosF = bpF;
                        long avgLastPosB = bpB;
                        
                        long chrPosStartF = (((long)chrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
                        long chrPosStopF = (((long)chrF<<28)+bpF)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
                        long chrPosStartB = (((long)chrB<<28)+bpB)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
                        long chrPosStopB = (((long)chrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

                        int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
                        int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
                        
                        String strAnnoF = "null";
                        String strAnnoB = "null";
                        if(annoGroupIndexF >= 0){
                            Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                            strAnnoF = annoF.toString();
                        }
                        if(annoGroupIndexB >= 0){
                            Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                            strAnnoB = annoB.toString();
                        }
                        
                        
                        /****************/
                        
                        writer.write("Group "+count);
                        writer.write("\tFront Break point = " + chrF +":"+bpF);
                        writer.write("\tAnnotation = " + strAnnoF);
                        writer.write("\tBack Break point = " + chrB +":"+bpB);
                        writer.write("\tAnnotation = " + strAnnoB);
    //                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            writer.write(var.virtualSequence());
                            writer.write("\n");
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
//                            writer.write(" || ");
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
//                            writer.write("\n");
                        }
                        count++;
                    }
                }  
            }
        }
        if(varType == 'I'){
            int count = 1;
            writer.write("Variation type : Indel\n");
            Set set = this.coverageMapIndel.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
                long bpFCode = (long)iterKey.next();
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                
                Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
                Map<Long,ArrayList<Integer>> sumNumMatchMapII = this.sumNumMatchCoverageMapFusion.get(bpFCode);
                Set setII = coverageMapII.keySet();
                Iterator iterKeyII = setII.iterator();
                while(iterKeyII.hasNext()){
                    long bpBCode = (long)iterKeyII.next();
                    long bpB = bpBCode&this.mask28bit;
                    int chrB = (int)(bpBCode>>28);
                    
                    /**
                     * Write Report Part
                     */
                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    ArrayList<Integer> sumNumMatchList = sumNumMatchMapII.get(bpBCode);
                    if(coverageList.size()>coverageThreshold){
                        /**
                         * pick one variation to get indel information 
                         */
                        Variation tempVar = coverageList.get(0);
                        String indelType = tempVar.getIndelType();
                        long indelBase = tempVar.getIndelBase();
                        /***/
                                               
                        /**
                         * Calculate average numMatch for Front and Back part and Map to annotation binary tree for gene information
                         * create code chrPosStart and chrPosStop both front and back (use for binary search)
                         */
                        
                        int avgNumMatchF = sumNumMatchList.get(0)/coverageList.size();
                        int avgNumMatchB = sumNumMatchList.get(1)/coverageList.size();
                        
//                        long avgIniPosF = bpF-avgNumMatchF;
//                        long avgLastPosB = bpB+avgNumMatchB;
                        
                        long avgIniPosF = bpF;
                        long avgLastPosB = bpB;
                        
                        long chrPosStartF = (((long)chrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
                        long chrPosStopF = (((long)chrF<<28)+bpF)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
                        long chrPosStartB = (((long)chrB<<28)+bpB)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
                        long chrPosStopB = (((long)chrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

                        int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
                        int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
                        
                        String strAnnoF = "null";
                        String strAnnoB = "null";
                        if(annoGroupIndexF >= 0){
                            Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                            strAnnoF = annoF.toString();
                        }
                        if(annoGroupIndexB >= 0){
                            Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                            strAnnoB = annoB.toString();
                        }
                        
                        /****************/
                        
                        writer.write("Group "+count);
                        writer.write("\tIndel Type : " + indelType);
                        writer.write("\tIndel Base : " + indelBase);
                        writer.write("\tFront Break point : " + chrF +":"+bpF);
                        writer.write("\tAnnotation = " + strAnnoF);
                        writer.write("\tBack Break point : " + chrB +":"+bpB);
                        writer.write("\tAnnotation = " + strAnnoB);
                        writer.write("\tCoverage : " + coverageList.size());
                        writer.write("\n");
                        for(int i=0;i<coverageList.size();i++){
                            Variation var = coverageList.get(i);
                            
                            writer.write(var.virtualSequence());
                            writer.write("\n");
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF,var.readLengthF));
//                            writer.write(" || ");
//                            writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB,var.readLengthB));
//                            writer.write("\n");
                        }
                        count++;
                    }
                }
                
            }
            
        }
        if(varType == 'S'){
            int count = 1;
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +","+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);
                        
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            
        }
        if(varType == 'O'){
            
        }

        writer.flush();
        writer.close();
    }
    
    public void writeVariantSortedCoverageVirtualizeReportToFile(String nameFile , int coverageThreshold , char varType) throws IOException{
        /**
        * Suitable for version 3 data structure (data structure that has iniIdx in its)
        * write result to file format for virtualize variant report in sorted order (high to low)
        * If varTyp is indel write into 3 separate file (small deletion, large deletion and insertion)
        * "Caution : this function must be called after sortCoverage[Indel or Fusion] function"
        */
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        String sdFile = "";
        String idFile = "";
        String ldFile = "";
        
        switch (varType) {
            case 'F':
                filename = dummy[0]+"_match"+this.percentMatch+"_SortedVirtualCoverageReport_Fusion"+".txt";                
                break;
            case 'I':
                sdFile = dummy[0]+"_match"+this.percentMatch+"_SortedVirtualCoverageReport_smallDelete"+".txt";
                idFile = dummy[0]+"_match"+this.percentMatch+"_SortedVirtualCoverageReport_smallInsert"+".txt";
                ldFile = dummy[0]+"_match"+this.percentMatch+"_SortedVirtualCoverageReport_largeDelete"+".txt";
                break;
            case 'S':
                filename = dummy[0]+"_match"+this.percentMatch+"_SortedVirtualCoverageReport_SNP"+".txt";                
                break;
            default:
                break;
        }

        PrintStream ps;
        FileWriter writer;
        FileWriter writerSD;    // small indel
        FileWriter writerID;    // small insertion
        FileWriter writerLD;    // large indel
        
        if(varType == 'F'){
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
            writer.write("Variation type : Fusion\n");
            for(int i=0;i<this.sortedCoverageArrayFusion.size();i++){
                if(this.sortedCoverageArrayFusion.size() == 0){
                    break;
                }
                ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayFusion.get(i);
                long bpFCode = dummyCoverageList.get(1);
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                long bpBCode = dummyCoverageList.get(2);
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                
                ArrayList<Variation> coverageList = getFusionCoverageList(bpFCode,bpBCode);
                
                if(coverageList.size()>coverageThreshold){
                    writer.write("Group "+count);
                    writer.write("\tFront Break point : " + chrF +":"+bpF);
                    writer.write("\tBack Break point : " + chrB +":"+bpB); 
//                    ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");
                    for(int j=0;j<coverageList.size();j++){
                        Variation var = coverageList.get(j);

                        writer.write(var.virtualSequence());
                        writer.write("\n");
                    }
                    count++;
                }
            }
            writer.flush();
            writer.close();
        }
        if(varType == 'I'){
            
            writerSD = new FileWriter(sdFile);  //small delete File writer
            writerID = new FileWriter(idFile);  //small insert File writer
            writerLD = new FileWriter(ldFile);  //large delete File writer
            
            int count = 1;
            int countSD = 1;
            int countID = 1;
            int countLD = 1;
//            writer.write("Variation type : Indel\n");
            
            for(int i=0;i<this.sortedCoverageArrayIndel.size();i++){
                if(this.sortedCoverageArrayIndel.size() == 0){
                    break;
                }
                ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayIndel.get(i);
                long bpFCode = dummyCoverageList.get(1);
                long bpF = bpFCode&this.mask28bit;
                int chrF = (int)(bpFCode>>28);
                long bpBCode = dummyCoverageList.get(2);
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = getIndelCoverageList(bpFCode,bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                     * pick one variation to get indel information 
                     */
                    Variation tempVar = coverageList.get(0);
                    String indelType = tempVar.getIndelType();
                    long indelBase = tempVar.getIndelBase();
                    /***/
                    if(indelType.equals("delete")){
                        writerSD.write("Group "+countSD);
                        writerSD.write("\tIndel Type : " + indelType);
                        writerSD.write("\tIndel Base : " + indelBase);
                        writerSD.write("\tFront Break point : " + chrF +":"+bpF);
                        writerSD.write("\tBack Break point : " + chrB +":"+bpB);
                        writerSD.write("\tCoverage : " + coverageList.size());
                        writerSD.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);

                            writerSD.write(var.virtualSequence());
                            writerSD.write("\n");
                        }
                        countSD++;
//                        File indelFile = new File(sdFile); //File object        
//                        if(indelFile.exists()){
//                            writer = new FileWriter(sdFile,true);
//                        }else{
//                            writer = new FileWriter(sdFile);
//                        }   
                    }else if(indelType.equals("insert")){
                        writerID.write("Group "+countID);
                        writerID.write("\tIndel Type : " + indelType);
                        writerID.write("\tIndel Base : " + indelBase);
                        writerID.write("\tFront Break point : " + chrF +":"+bpF);
                        writerID.write("\tBack Break point : " + chrB +":"+bpB);
                        writerID.write("\tCoverage : " + coverageList.size());
                        writerID.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);

                            writerID.write(var.virtualSequence());
                            writerID.write("\n");
                        }
                        countID++;
                        
                    }else{
                        // large indel
                        writerLD.write("Group "+countLD);
                        writerLD.write("\tIndel Type : " + indelType);
                        writerLD.write("\tIndel Base : " + indelBase);
                        writerLD.write("\tFront Break point : " + chrF +":"+bpF);
                        writerLD.write("\tBack Break point : " + chrB +":"+bpB);
                        writerLD.write("\tCoverage : " + coverageList.size());
                        writerLD.write("\n");
                        for(int j=0;j<coverageList.size();j++){

                            Variation var = coverageList.get(j);

                            writerLD.write(var.virtualSequence());
                            writerLD.write("\n");
                        }
                        countLD++;
                    }
                }  
            }
            writerSD.flush();
            writerSD.close();
            
            writerID.flush();
            writerID.close();
            
            writerLD.flush();
            writerLD.close();      
        }
        if(varType == 'S'){
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
            writer.write("Variation type : SNP and other miss match\n");
            Set set = this.coverageMapSNP.keySet();
            Iterator iterKey = set.iterator();
            while(iterKey.hasNext()){
//                long bpFCode = (long)iterKey.next();
//                long bpF = bpFCode&this.mask28bit;
//                int chrF = (int)bpFCode>>28;
                long bpBCode = (long)iterKey.next();
                long bpB = bpBCode&this.mask28bit;
                int chrB = (int)(bpBCode>>28);
                
                ArrayList<Variation> coverageList = this.coverageMapSNP.get(bpBCode);
                if(coverageList.size()>coverageThreshold){
                    /**
                    * Write Report Part
                    */
                    writer.write("Group "+count);
                    writer.write("\tBack Break point : " + chrB +":"+bpB);  
                    writer.write("\tCoverage : " + coverageList.size());
                    writer.write("\n");

                    for(int i=0;i<coverageList.size();i++){
                        /**
                         * Write Report Part
                         */
                        Variation var = coverageList.get(i);

                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrF,var.iniPosF,var.lastPosF,var.greenF,var.yellowF,var.orangeF,var.redF,var.strandF,var.iniIndexF,var.readNameF,var.snpFlagF,var.iniBackFlagF));
                        writer.write(" || ");
                        writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%d,%d", var.numChrB,var.iniPosB,var.lastPosB,var.greenB,var.yellowB,var.orangeB,var.redB,var.strandB,var.iniIndexB,var.readNameB,var.snpFlagB,var.iniBackFlagB));
                        writer.write("\n");
                    }
                    count++;
                }
            }
            writer.flush();
            writer.close();
            
        }
        if(varType == 'O'){
            
        }

        

    }
    
    public static boolean containsVariation(ArrayList<Variation> list, Variation inVar) {
        /**
         * Utility function for checking the inVar is contain in ArrayList of variation or not
         * Check by read name
         */
        for (Variation object : list) {
            if (object.readNameF.equals(inVar.readNameF)) {
                return true;
            }
        }
        return false;
    }
    
    public void sortCoverageIndel(){
        /**
         * This function will create a Array of coverageIndel that sort by high to low (descending)
         * This array will be use as index to call coverage information from map in descending order
         */
        
        Set set = this.coverageMapIndel.keySet();
        Iterator iterKey = set.iterator();
        while(iterKey.hasNext()){
            long bpFCode = (long)iterKey.next();
//            long bpF = bpFCode&this.mask28bit;
//            int chrF = (int)(bpFCode>>28);

            Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapIndel.get(bpFCode);
            Set setII = coverageMapII.keySet();
            Iterator iterKeyII = setII.iterator();
            while(iterKeyII.hasNext()){
                long bpBCode = (long)iterKeyII.next();               
//                long bpB = bpBCode&this.mask28bit;
//                int chrB = (int)(bpBCode>>28);

                ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                long numCoverage = coverageList.size();
                ArrayList<Long> signatureArray = new ArrayList();
                signatureArray.add(numCoverage);
                signatureArray.add(bpFCode);
                signatureArray.add(bpBCode);
                
                this.sortedCoverageArrayIndel.add(signatureArray);
            }    
        }
        
        // Begin sort Array by implement comparator to force it to sort from element 0 of array which is a numCoverage
        Collections.sort(this.sortedCoverageArrayIndel, new Comparator<ArrayList<Long>>() {
            @Override
            public int compare(final ArrayList<Long> a1, final ArrayList<Long> a2) {
                final long sigArr1 = a1.get(0);
                final long sigArr2 = a2.get(0);
                return (int)(sigArr2-sigArr1); // sort by descending (from high to low)   if yo want to sort by ascending just switch btw 1 and 2
            }
        });
        /**********************************************************************/   
    }
    
    public void sortCoverageFusion(){
        /**
         * This function will create a Array of coverageFusion that sort by high to low (descending)
         * This array will be use as index to call coverage information from map in descending order
         */
        
        Set set = this.coverageMapFusion.keySet();
        Iterator iterKey = set.iterator();
        while(iterKey.hasNext()){
            long bpFCode = (long)iterKey.next();
//            long bpF = bpFCode&this.mask28bit;
//            int chrF = (int)(bpFCode>>28);

            Map<Long,ArrayList<Variation>> coverageMapII = this.coverageMapFusion.get(bpFCode);
            Set setII = coverageMapII.keySet();
            Iterator iterKeyII = setII.iterator();
            while(iterKeyII.hasNext()){
                long bpBCode = (long)iterKeyII.next();               
//                long bpB = bpBCode&this.mask28bit;
//                int chrB = (int)(bpBCode>>28);

                ArrayList<Variation> coverageList = coverageMapII.get(bpBCode);
                long numCoverage = coverageList.size();
                ArrayList<Long> signatureArray = new ArrayList();
                signatureArray.add(numCoverage);
                signatureArray.add(bpFCode);
                signatureArray.add(bpBCode);
                
                this.sortedCoverageArrayFusion.add(signatureArray);
            }    
        }
        
        // Begin sort Array by implement comparator to force it to sort from element 0 of array which is a numCoverage
        Collections.sort(this.sortedCoverageArrayFusion, new Comparator<ArrayList<Long>>() {
            @Override
            public int compare(final ArrayList<Long> a1, final ArrayList<Long> a2) {
                final long sigArr1 = a1.get(0);
                final long sigArr2 = a2.get(0);
                return (int)(sigArr2-sigArr1); // sort by descending (from high to low)   if yo want to sort by ascending just switch btw 1 and 2
            }
        });
        /**********************************************************************/
    }
    
    public ArrayList<Variation> getIndelCoverageList(long signature1, long signature2){
        Map<Long,ArrayList<Variation>> mapII = this.coverageMapIndel.get(signature1);
        ArrayList<Variation> indelCoverageList = mapII.get(signature2);
        return indelCoverageList;
    }
    
    public ArrayList<Variation> getFusionCoverageList(long signature1, long signature2){
        Map<Long,ArrayList<Variation>> mapII = this.coverageMapFusion.get(signature1);
        ArrayList<Variation> fusionCoverageList = mapII.get(signature2);
        return fusionCoverageList;
    }
    
    public void exportNewIndelEventToBed12(String nameFile, String inIndelType, int minPickCoverage) throws IOException{
        /**
         * This function will pick top10/20 or user define of specific indelType (SI = small Insertion, SD = small deletion, LI = large Indel)
         * Calculate range of new reference from breakpoint with user define length EX bpF=100 and bpB=200 want to create reference 200 bp length 
         * We use bed12 format to represent the position that we want to extract from reference. Bed12 format of this example should be
         *  chr1 0 300 readname 0 + 0 300 0 2 100,100 0,200
         * 
         * The position will be extend from breakpoint front for 300 base on the left and breakpoint back for 300 base on the right
         * 
         * And export to File in format bed for create new reference with bedtools
         * 
         * ##Not finish
         */
        String[] dummy = nameFile.split("\\.");
        String filename = "";
        FileWriter writer;        
        int extendSize = 300;
        
        String indelType = "";
        if(inIndelType.equals("SI")){
            indelType = "insert";
            filename = dummy[0]+"newSmallInsert"+".bed";
        }else if(inIndelType.equals("SD")){
            indelType = "delete";
            filename = dummy[0]+"_newSmallDelete"+".bed";
        }else if(inIndelType.equals("LI")){
            indelType = "large indel";
            filename = dummy[0]+"_newLargeDelete"+".bed";
        }
        
        writer = new FileWriter(filename);  
              
        for(int i=0;i<this.sortedCoverageArrayIndel.size();i++){
            if(this.sortedCoverageArrayIndel.size() == 0){
                break;
            }

            ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayIndel.get(i);
            long bpFCode = dummyCoverageList.get(1);
            long bpF = bpFCode&this.mask28bit;
            int chrF = (int)(bpFCode>>28);
            long bpBCode = dummyCoverageList.get(2);
            long bpB = bpBCode&this.mask28bit;
            int chrB = (int)(bpBCode>>28);

            ArrayList<Variation> coverageList = getIndelCoverageList(bpFCode,bpBCode);
            
            if(coverageList.size() >= minPickCoverage){
                Variation dummyVariation = coverageList.get(0);
                if(dummyVariation.getIndelType().equals(indelType)){                       
                    String bed12Format = dummyVariation.exportBed12Ref(extendSize);
                    writer.write(bed12Format);
                    writer.write("\n");
                } 
            }else{
                break;
            }   
        }
        
        writer.flush();
        writer.close();
    }
    
    public void createReferenceFromNovelIndelResult(String nameFile, String refFile, String refIdxFile, String readFile, String readIdxFile, String inIndelType, int minPickCoverage, int inExtendSize) throws IOException{
        /**
         * This function will pick group that past user define coverage of specific indelType (SI = small Insertion, SD = small deletion, LI = large Indel)
         * It will pick first variation pattern of each group to create new reference 
         * Extending left and right of the picked read sequence by cut the left wing and right wing sequence from reference
         * Concatenate left read seq and right wing together
         * Export into new reference fasta file
         * 
         * EX read has 100 base long and extendSize is 200 
         * at the end you will got new ref with 200 + 100 + 200 = 500 base long
         * 
         * ##Not finish
         */
        RandomAccessFile rbRef = new RandomAccessFile(refFile,"r");
        RandomAccessFile rbSample = new RandomAccessFile(readFile,"r");
        
        RefFaIndex refFaIdx = new RefFaIndex(refIdxFile);
        SampleFaIndex sampleFaIdx = new SampleFaIndex(readIdxFile);
//        /**
//         * read reference fa index file
//         * byte will count in zero based position (first index start is zero)
//         */
//        RandomAccessFile rbRefIdx = new RandomAccessFile(refIdxFile,"r");
//        String line;
//        String name = "";
//        long totalLen = 0;
//        long offset = 0;        // offset or in another mean is pointer that point to the number of byte of first base of sequence
//        long lineBase = 0;      // number of base on each line (some sequence may cover many of line ex 200 base long may cover 10 line if it wirte 20 base per line)
//        long lineWidth = 0;     // number of byte in each line (include the new line) EX 1 byte per base 20 base per line so lineWidth = 21 byte (20 base + newline)
//        Map<String,ArrayList<Long>> refFaIdx = new LinkedHashMap();     // map that store chromosome name as key and ArrayList of index info as value
//        
//        while ((line = rbRefIdx.readLine()) != null) {
//            String[] linePortion = line.split("\t");
//            name = linePortion[0];
//            totalLen = Long.parseLong(linePortion[1]);
//            offset = Long.parseLong(linePortion[2]);
//            lineBase = Long.parseLong(linePortion[3]);
//            lineWidth = Long.parseLong(linePortion[4]);
//            ArrayList<Long> listIdxInfo = new ArrayList();
//            listIdxInfo.add(totalLen);
//            listIdxInfo.add(offset);
//            listIdxInfo.add(lineBase);
//            listIdxInfo.add(lineWidth);
//            
//            refFaIdx.put(name, listIdxInfo);
//        }
//  
//        /******************************/
//        
//        /**
//         * read read sample fa index file
//         * 
//         */
//        RandomAccessFile rbSampleIdx = new RandomAccessFile(refIdxFile,"r");
//        Map<String,ArrayList<Long>> sampleFaIdx = new LinkedHashMap();     // map that store chromosome name as key and ArrayList of index info as value
//        
//        while ((line = rbSampleIdx.readLine()) != null) {
//            String[] linePortion = line.split("\t");
//            name = linePortion[0];
//            totalLen = Long.parseLong(linePortion[1]);
//            offset = Long.parseLong(linePortion[2]);
//            lineBase = Long.parseLong(linePortion[3]);
//            lineWidth = Long.parseLong(linePortion[4]);
//            ArrayList<Long> listIdxInfo = new ArrayList();
//            listIdxInfo.add(totalLen);
//            listIdxInfo.add(offset);
//            listIdxInfo.add(lineBase);
//            listIdxInfo.add(lineWidth);
//            
//            sampleFaIdx.put(name, listIdxInfo);
//        }
//        
//        /*******************************/
        
        String[] dummy = readFile.split("\\.");
        Path path = Paths.get(readFile);
        String fileName = path.getFileName().toString();
        String[] dummy2 = fileName.split("_");
        String sampleName = dummy2[0];
        String filename = "";
        FileWriter writer;        
        int extendSize = inExtendSize;
        int groupCount = 0;
        
        String indelType = "";
        if(inIndelType.equals("SI")){
            indelType = "insert";
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceSmallInsert"+".fa";
        }else if(inIndelType.equals("SD")){
            indelType = "delete";
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceSmallDelete"+".fa";
        }else if(inIndelType.equals("LI")){
            indelType = "large indel";
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceLargeDelete"+".fa";
        }
        
        writer = new FileWriter(filename);  
        
        /**
         * Start loop in array of sort coverage list
         */
        for(int i=0;i<this.sortedCoverageArrayIndel.size();i++){
            if(this.sortedCoverageArrayIndel.size() == 0){
                break;
            }

            ArrayList<Long> dummyCoverageList = this.sortedCoverageArrayIndel.get(i);
            long bpFCode = dummyCoverageList.get(1);
            long bpF = bpFCode&this.mask28bit;
            int chrF = (int)(bpFCode>>28);
            long bpBCode = dummyCoverageList.get(2);
            long bpB = bpBCode&this.mask28bit;
            int chrB = (int)(bpBCode>>28);

            ArrayList<Variation> coverageList = getIndelCoverageList(bpFCode,bpBCode);      // get coverage list form signature variabl (bpFCode and bpBCode)
            
            if(coverageList.size() >= minPickCoverage){
                /**
                 * pick first variation object to check for indel type and get 
                 */
                Variation dummyVariation = coverageList.get(0);
                if(dummyVariation.getIndelType().equals(indelType)){
                    groupCount++;

                    String leftWingSeq = "";
                    String rightWingSeq = "";
                    String readSeq = "";
                    String newRefSeq = "";
                    String strand = dummyVariation.strandF + dummyVariation.strandB;
                    
                    long[] signatureValue =  dummyVariation.getSignatureForCreateRef(extendSize);
                    long startLeftWing = signatureValue[0];
                    long startRightWing = signatureValue[1];
                    int iniIndexMatch = (int)signatureValue[2];
                    int lastIndexMatch = (int)signatureValue[3];
                    int unmatchFront = (int)signatureValue[4];
                    int unmatchBack = (int)signatureValue[5];
                    int numChrF = dummyVariation.numChrF;
                    int numChrB = dummyVariation.numChrB;
                    
                    String chrNameF = "chr"+numChrF;
                    String chrNameB = "chr"+numChrB;
                    
                    long compensateLeftLineByte = startLeftWing/refFaIdx.getLineBase(chrNameF);    // number of compensate byte that has to be adding back to the pointer (cause from newline byte in file)
                    long compensateRightLineByte = startRightWing/refFaIdx.getLineBase(chrNameB);  // number of compensate byte that has to be adding back to the pointer (cause from newline byte in file)
                    
                    long iniLeftWingPointer = ((startLeftWing + refFaIdx.getOffSet(chrNameF)))+compensateLeftLineByte;
                    long iniRightWingPointer = ((startRightWing + refFaIdx.getOffSet(chrNameB)))+compensateRightLineByte;
                    
                    long numLineRead = (extendSize/refFaIdx.getLineBase(chrNameF))+2;
                    
                    long readPointer = sampleFaIdx.getOffSet(dummyVariation.readNameF);
                    
                    rbRef.seek(iniLeftWingPointer);
                    String dummyLeftWingSeq = "";
                    for(int num = 0;num<numLineRead;num++){
                        dummyLeftWingSeq = dummyLeftWingSeq + rbRef.readLine();
                    }
                    leftWingSeq = dummyLeftWingSeq.substring(0, (extendSize+unmatchFront));
                    
                    rbRef.seek(iniRightWingPointer);
                    String dummyRightWingSeq = "";
                    for(int num = 0;num<numLineRead;num++){
                        dummyRightWingSeq = dummyRightWingSeq + rbRef.readLine();
                    }
                    rightWingSeq = dummyRightWingSeq.substring(0,(extendSize+unmatchBack));
                    
                    rbSample.seek(readPointer);
                    readSeq = rbSample.readLine().substring(iniIndexMatch, lastIndexMatch+1);
                    
                    if(strand.equals("++")){
                        newRefSeq = leftWingSeq + readSeq + rightWingSeq;
                    }else if(strand.equals("--")){
                        String leftWingInvSeq = SequenceUtil.inverseSequence(leftWingSeq);                                // Do invert sequence (ATCG => GCTA)
                        String leftWingCompSeq = SequenceUtil.createComplimentV2(leftWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                        
                        String rightWingInvSeq = SequenceUtil.inverseSequence(rightWingSeq);                                // Do invert sequence (ATCG => GCTA)
                        String rightWingCompSeq = SequenceUtil.createComplimentV2(rightWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                        
                        newRefSeq = leftWingCompSeq + readSeq + rightWingCompSeq;
                        
                    }else if(strand.equals("+-")){
                        String rightWingInvSeq = SequenceUtil.inverseSequence(rightWingSeq);                                // Do invert sequence (ATCG => GCTA)
                        String rightWingCompSeq = SequenceUtil.createComplimentV2(rightWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                        
                        newRefSeq = leftWingSeq + readSeq + rightWingCompSeq;
                    }else if(strand.equals("-+")){
                        String leftWingInvSeq = SequenceUtil.inverseSequence(leftWingSeq);                                // Do invert sequence (ATCG => GCTA)
                        String leftWingCompSeq = SequenceUtil.createComplimentV2(leftWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                        
                        newRefSeq = leftWingCompSeq + readSeq + rightWingSeq;
                    }
                    
                    String nameRef = "group"+groupCount+"_"+dummyVariation.readNameF;
                    
                    writer.write(">"+nameRef);
                    writer.write("\n");
                    writer.write(newRefSeq);
                    writer.write("\n");
                } 
            }else{
                break;
            }   
        }
        
        writer.flush();
        writer.close();
        
    }
    
    public void createReferenceFromNovelIndelResult_VariationV2(String readFile, String refFile, String refIdxFile, String inIndelType, int minPickCoverage, int inExtendSize) throws IOException{
        /**
         * This function will pick group that past user define coverage of specific indelType (SI = small Insertion, SD = small deletion, LI = large Indel)
         * It will pick first variation pattern of each group to create new reference 
         * Extending left and right of the picked read sequence by cut the left wing and right wing sequence from reference
         * Concatenate left read seq and right wing together
         * Export into new reference fasta file
         * 
         * EX read has 100 base long and extendSize is 200 
         * at the end you will got new ref with 200 + 100 + 200 = 500 base long
         * 
         * ##Not finish
         */
        RandomAccessFile rbRef = new RandomAccessFile(refFile,"r");
        
        RefFaIndex refFaIdx = new RefFaIndex(refIdxFile);

//        /**
//         * read reference fa index file
//         * byte will count in zero based position (first index start is zero)
//         */
//        RandomAccessFile rbRefIdx = new RandomAccessFile(refIdxFile,"r");
//        String line;
//        String name = "";
//        long totalLen = 0;
//        long offset = 0;        // offset or in another mean is pointer that point to the number of byte of first base of sequence
//        long lineBase = 0;      // number of base on each line (some sequence may cover many of line ex 200 base long may cover 10 line if it wirte 20 base per line)
//        long lineWidth = 0;     // number of byte in each line (include the new line) EX 1 byte per base 20 base per line so lineWidth = 21 byte (20 base + newline)
//        Map<String,ArrayList<Long>> refFaIdx = new LinkedHashMap();     // map that store chromosome name as key and ArrayList of index info as value
//        
//        while ((line = rbRefIdx.readLine()) != null) {
//            String[] linePortion = line.split("\t");
//            name = linePortion[0];
//            totalLen = Long.parseLong(linePortion[1]);
//            offset = Long.parseLong(linePortion[2]);
//            lineBase = Long.parseLong(linePortion[3]);
//            lineWidth = Long.parseLong(linePortion[4]);
//            ArrayList<Long> listIdxInfo = new ArrayList();
//            listIdxInfo.add(totalLen);
//            listIdxInfo.add(offset);
//            listIdxInfo.add(lineBase);
//            listIdxInfo.add(lineWidth);
//            
//            refFaIdx.put(name, listIdxInfo);
//        }
//  
//        /******************************/
//        
//        /**
//         * read read sample fa index file
//         * 
//         */
//        RandomAccessFile rbSampleIdx = new RandomAccessFile(refIdxFile,"r");
//        Map<String,ArrayList<Long>> sampleFaIdx = new LinkedHashMap();     // map that store chromosome name as key and ArrayList of index info as value
//        
//        while ((line = rbSampleIdx.readLine()) != null) {
//            String[] linePortion = line.split("\t");
//            name = linePortion[0];
//            totalLen = Long.parseLong(linePortion[1]);
//            offset = Long.parseLong(linePortion[2]);
//            lineBase = Long.parseLong(linePortion[3]);
//            lineWidth = Long.parseLong(linePortion[4]);
//            ArrayList<Long> listIdxInfo = new ArrayList();
//            listIdxInfo.add(totalLen);
//            listIdxInfo.add(offset);
//            listIdxInfo.add(lineBase);
//            listIdxInfo.add(lineWidth);
//            
//            sampleFaIdx.put(name, listIdxInfo);
//        }
//        
//        /*******************************/
        
//        String[] dummy = readFile.split("\\.");
        Path path = Paths.get(readFile);
        String fileName = path.getFileName().toString();
        String[] dummy2 = fileName.split("_");
        String sampleName = dummy2[0];
        String filename = "";
        FileWriter writer;        
        int extendSize = inExtendSize;
        int groupCount = 0;
        ArrayList<SVGroup> variationList = new ArrayList();
        
        String indelType = "";
        if(inIndelType.equals("TD")){
            variationList = this.tandemList;
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceTandem"+".fa";
        }else if(inIndelType.equals("ID")){
            variationList = this.indelList;
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceIndel"+".fa";
        }else if(inIndelType.equals("IC")){
            variationList = this.intraTransList;
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceIntraTrans"+".fa";
        }else if(inIndelType.equals("IT")){
            variationList = this.interTransList;
            filename = path.getParent().toString()+File.separator+sampleName+"_newReferenceInterTrans"+".fa";
        }
        
        writer = new FileWriter(filename);  
        
        /**
         * Start loop in array of sort coverage list
         */
        for(int i=0;i<variationList.size();i++){
            if(variationList.size() == 0){
                break;
            }
            
            SVGroup dummySVGroup = variationList.get(i);
            
            if(dummySVGroup.getNumCoverage() >= minPickCoverage){
                /**
                 * pick first variation object to check for indel type and get 
                 */
                VariationV2 dummyVariationV2 = dummySVGroup.getVarList().get(0);        // pick first read of this group
                
//                Variation dummyVariation = coverageList.get(0);
                
                groupCount++;

                String leftWingSeq = "";
                String rightWingSeq = "";
                String readSeq = "";
                String newRefSeq = "";
                String strand = dummyVariationV2.strandF+ "" + dummyVariationV2.strandB;

                long[] signatureValue =  dummyVariationV2.getSignatureForCreateRef(extendSize);
                long startLeftWing = signatureValue[0];
                long startRightWing = signatureValue[1];
                int iniIndexMatch = (int)signatureValue[2];
                int lastIndexMatch = (int)signatureValue[3];
                int unmatchFront = (int)signatureValue[4];
                int unmatchBack = (int)signatureValue[5];

                String chrNameF = dummyVariationV2.chrF;
                String chrNameB = dummyVariationV2.chrB;

                long compensateLeftLineByte = startLeftWing/refFaIdx.getLineBase(chrNameF);    // number of compensate byte that has to be adding back to the pointer (cause from newline byte in file)
                long compensateRightLineByte = startRightWing/refFaIdx.getLineBase(chrNameB);  // number of compensate byte that has to be adding back to the pointer (cause from newline byte in file)

                long iniLeftWingPointer = ((startLeftWing + refFaIdx.getOffSet(chrNameF)))+compensateLeftLineByte;
                long iniRightWingPointer = ((startRightWing + refFaIdx.getOffSet(chrNameB)))+compensateRightLineByte;

                long numLineRead = (extendSize/refFaIdx.getLineBase(chrNameF))*2;

//                long readPointer = sampleFaIdx.getOffSet(dummyVariation.readNameF);

                rbRef.seek(iniLeftWingPointer);
                String dummyLeftWingSeq = "";
                for(int num = 0;num<numLineRead;num++){
                    dummyLeftWingSeq = dummyLeftWingSeq + rbRef.readLine();
                }
//                if(dummyVariationV2.readID == 71575276){
//                    System.out.println(dummyVariationV2.readID);
//                    System.out.println("long="+dummyLeftWingSeq.length());
//                    System.out.println(extendSize+unmatchFront);
//                }
                
                leftWingSeq = dummyLeftWingSeq.substring(0, (extendSize+unmatchFront));

                rbRef.seek(iniRightWingPointer);
                String dummyRightWingSeq = "";
                for(int num = 0;num<numLineRead;num++){
                    dummyRightWingSeq = dummyRightWingSeq + rbRef.readLine();
                }
                rightWingSeq = dummyRightWingSeq.substring(0,(extendSize+unmatchBack));

//                rbSample.seek(readPointer);
//                readSeq = rbSample.readLine().substring(iniIndexMatch, lastIndexMatch+1);
                readSeq = dummyVariationV2.getReadSeq().substring(iniIndexMatch, lastIndexMatch+1);
                
                if(strand.equals("00")){
                    newRefSeq = leftWingSeq + readSeq + rightWingSeq;
                }else if(strand.equals("11")){
                    String leftWingInvSeq = SequenceUtil.inverseSequence(leftWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String leftWingCompSeq = SequenceUtil.createComplimentV2(leftWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    String rightWingInvSeq = SequenceUtil.inverseSequence(rightWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String rightWingCompSeq = SequenceUtil.createComplimentV2(rightWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    newRefSeq = leftWingCompSeq + readSeq + rightWingCompSeq;

                }else if(strand.equals("01")){
                    String rightWingInvSeq = SequenceUtil.inverseSequence(rightWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String rightWingCompSeq = SequenceUtil.createComplimentV2(rightWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    newRefSeq = leftWingSeq + readSeq + rightWingCompSeq;
                }else if(strand.equals("10")){
                    String leftWingInvSeq = SequenceUtil.inverseSequence(leftWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String leftWingCompSeq = SequenceUtil.createComplimentV2(leftWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    newRefSeq = leftWingCompSeq + readSeq + rightWingSeq;
                }

                String nameRef = "group"+groupCount+"_"+dummyVariationV2.readID;

                writer.write(">"+nameRef);
                writer.write("\n");
                writer.write(newRefSeq);
                writer.write("\n");   
            }else{
                break;
            }   
        }
        
        writer.flush();
        writer.close();
        
    }
    
    
    public void classifyRoughSVType(){
        /**
         * classify sv type of each sv candidate group
         * We hasve 4 rough type of SV
         *  1. tandem
         *  2. indel
         *  3. intraTrans
         *  4. interTrans
         * 
         * (extra** In order to get true type of SV, all 4 type hasve to be pass through relation analysis algorithm (other function))
         */
//        int count=0;
        for(Map.Entry<String, Map<String,SVGroup>> entry1 : this.coverageMap.entrySet()){
            Map<String,SVGroup> coverageMapII = entry1.getValue();
            
            for(Map.Entry<String,SVGroup> entry2 : coverageMapII.entrySet()){
                SVGroup svGroup = entry2.getValue();
//                System.out.println(count++);
//                if(count==26){
//                    System.out.println();
//                }
                if(svGroup.getSVType().equals("tandem")){
                    this.tandemList.add(svGroup);
                }else if(svGroup.getSVType().equals("indel")){
                    this.indelList.add(svGroup);
                }else if(svGroup.getSVType().equals("intraTrans")){
                    this.intraTransList.add(svGroup);
                }else if(svGroup.getSVType().equals("interTrans")){
                    this.interTransList.add(svGroup);
                }
            }
        }
        
        // sort list of SVGroup by number of coverage fron high to low (descending) comparable has benn inplement in SVGroup object
        Collections.sort(this.tandemList);
        Collections.sort(this.indelList);
        Collections.sort(this.intraTransList);
        Collections.sort(this.interTransList);        
    }
    
    public void writeStructureVariantV2SortedCoverageReportToFile(String nameFile , int coverageThreshold) throws IOException{
        /**
        * Suitable for SVGroup object that contain VariaitonV2 object
        * write result to file format for variant report in sorted order (high to low)
        * 
        * "Caution : this function must be called after classifyRoughVariantReport function"
        */
        String tandemFile = nameFile+".tandem_cov.out";
        String indelFile = nameFile+".indel_cov.out";
        String intraTransFile = nameFile+".intraTrans_cov.out";
        String interTransFile = nameFile+".interTrans_cov.out";
        String groupCoverageReport = nameFile+".mrkDup.cov.out";
        /**
         * Write tandem cov report
         */
        FileWriter writer;
        File f = new File(groupCoverageReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(groupCoverageReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(groupCoverageReport);
        }
        int count = 1;
        for(int i=0;i<this.tandemList.size();i++){
            SVGroup svGroup = this.tandemList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\n");
            for(int j=0;j<varList.size();j++){
                VariationV2 dummyVar = varList.get(j);
                writer.write(dummyVar.toString());
                writer.write("\n");
            }
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write indel cov report
         */
//        f = new File(indelFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(indelFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(indelFile);
//        }
//        count = 1;
        for(int i=0;i<this.indelList.size();i++){
            SVGroup svGroup = this.indelList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\n");
            for(int j=0;j<varList.size();j++){
                VariationV2 dummyVar = varList.get(j);
                writer.write(dummyVar.toString());
                writer.write("\n");
            }
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write intraTrans cov report
         */
//        f = new File(intraTransFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(intraTransFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(intraTransFile);
//        }
//        count = 1;
        for(int i=0;i<this.intraTransList.size();i++){
            SVGroup svGroup = this.intraTransList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\n");
            for(int j=0;j<varList.size();j++){
                VariationV2 dummyVar = varList.get(j);
                writer.write(dummyVar.toString());
                writer.write("\n");
            }
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write interTrans cov report
         */
//        f = new File(interTransFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(interTransFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(interTransFile);
//        }
//        count = 1;
        for(int i=0;i<this.interTransList.size();i++){
            SVGroup svGroup = this.interTransList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\n");
            for(int j=0;j<varList.size();j++){
                VariationV2 dummyVar = varList.get(j);
                writer.write(dummyVar.toString());
                writer.write("\n");
            }
            count++;
        }
        
        writer.flush();
        writer.close();
        /***************************************************/
    }
    
    public void writeStructureVariantV2SortedCoverageGroupInfoReportToFile(String nameFile , int coverageThreshold) throws IOException{
        /**
        * Suitable for SVGroup object that contain VariaitonV2 object
        * write result to file format for variant report in sorted order (high to low)
        * 
        * "Caution : this function must be called after classifyRoughVariantReport function"
        */
        String groupCoverageReport = nameFile+".mrkDup.cov.ginfo.out";
        /**
         * Write tandem cov report
         */
        FileWriter writer;
        File f = new File(groupCoverageReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(groupCoverageReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(groupCoverageReport);
        }
        int count = 1;
        for(int i=0;i<this.tandemList.size();i++){
            SVGroup svGroup = this.tandemList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+"\t"+svGroup.toString());
            writer.write("\n");
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write indel cov report
         */
//        f = new File(indelFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(indelFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(indelFile);
//        }
//        count = 1;
        for(int i=0;i<this.indelList.size();i++){
            SVGroup svGroup = this.indelList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+"\t"+svGroup.toString());
            writer.write("\n");
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write intraTrans cov report
         */
//        f = new File(intraTransFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(intraTransFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(intraTransFile);
//        }
//        count = 1;
        for(int i=0;i<this.intraTransList.size();i++){
            SVGroup svGroup = this.intraTransList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+"\t"+svGroup.toString());
            writer.write("\n");
            count++;
        }
        
//        writer.flush();
//        writer.close();
        /***************************************************/
        /**
         * Write interTrans cov report
         */
//        f = new File(interTransFile); //File object        
//        if(f.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(interTransFile,true);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(interTransFile);
//        }
//        count = 1;
        for(int i=0;i<this.interTransList.size();i++){
            SVGroup svGroup = this.interTransList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+"\t"+svGroup.toString());
            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /***************************************************/
    }
}
