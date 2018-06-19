/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
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
    private ArrayList<SVGroup> deletionList;
    private ArrayList<SVGroup> intraTransList;
    private ArrayList<SVGroup> interTransList;
    private ArrayList<SVGroupPair> intraInsertionList;
    private ArrayList<SVGroupPair> intraInsertionList_outFilter;
    private ArrayList<SVGroupPair> interInsertionList;
    private ArrayList<SVGroup> sameChrSVGroup;
    private ArrayList<SVGroup> diffChrSVGroup; 
    private Map<ArrayList<SVGroup>,Boolean> intraTransPairList;
    private Map<ArrayList<SVGroup>,Boolean> interTransPairList;
    private ArrayList<SVGroup> chimericList;
    
    private Map<String,Long> refIndex;                              // store chrmosome name as key and assigned chr number in value (this map is very important in case of chr name is not a number we will asigned the number to it and store in this map for further reference)
    
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
        this.deletionList = new ArrayList();
        this.intraTransList = new ArrayList();
        this.interTransList = new ArrayList();
        this.interInsertionList = new ArrayList();
        this.intraInsertionList = new ArrayList();
        this.refIndex = new LinkedHashMap();
        this.sameChrSVGroup = new ArrayList();
        this.diffChrSVGroup = new ArrayList();
        this.intraTransPairList = new LinkedHashMap();
        this.interTransPairList = new LinkedHashMap();
        this.chimericList = new ArrayList();
        this.intraInsertionList_outFilter = new ArrayList();
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
         * This fuction use to group same SV event together
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
                    svGroup.addRefIndex(this.refIndex);
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
                svGroup.addRefIndex(this.refIndex);
            }            
        }
//        System.out.println("");
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
    
    public void sortInsertionByInsertSize(){
        Collections.sort(this.intraInsertionList,SVGroupPair.InsertSizeComparatorLowToHigh);
        Collections.sort(this.interInsertionList,SVGroupPair.InsertSizeComparatorLowToHigh);
    }
    
    public void sortInsertionByInsertJunctiontLowtoHigh(){
        Collections.sort(this.intraInsertionList,SVGroupPair.InsertJunctionComparatorLowToHigh);
        Collections.sort(this.interInsertionList,SVGroupPair.InsertJunctionComparatorLowToHigh);
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
        rbRef.close();
        rbSample.close();
        
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
         * Implement save single junction file
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
        String[] dummy2 = fileName.split("_unmap");
        String sampleName = dummy2[0];
        String filename = "";
        FileWriter writer;
        FileWriter writer2;
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
        rbRef.close();
        
    }
    
    public void createReferenceFromNovelIndelResult_VariationV3(String readFile, String refFile, String refIdxFile, String inIndelType, int minPickCoverage, int inExtendSize) throws IOException{
        /**
         * This function will pick group that past user define coverage of specific indelType (SI = small Insertion, SD = small deletion, LI = large Indel)
         * It will pick first variation pattern of each group to create new reference 
         * Extending left and right of the picked read sequence by cut the left wing and right wing sequence from reference
         * Concatenate left read seq and right wing together
         * Export into new reference fasta file (1 junction per file)
         * 
         * EX read has 100 base long and extendSize is 200 
         * at the end you will got new ref with 200 + 100 + 200 = 500 base long
         * 
         * Add ArrayList of string, the string represent the region on reference that we consider to create new reference. This String will write in to txt file to use as a region specific filter when filter bam file
         * This String will have format chromosomeName:leftPos-rightPos and have space delimit for multiple region
         * The leftPos and rightPos will calculate from One wing (left or right wing). the leftPos calculate from left/right wing - (extendSize*2) . the rightPos calculate from left/right wing + (extendSize*2) 
         * 
         * Implement save separate junction file
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
        String[] dummy2 = fileName.split("_unmap");
        String sampleName = dummy2[0];
        String filename = "";
        String regionFilterFileName = "";
        String regionFilter="";
        FileWriter writer;
        FileWriter writer2;
        int extendSize = inExtendSize;
        int groupCount = 0;
        long regionLeftWingLeft = 0;
        long regionLeftWingRight = 0;
        long regionRightWingLeft = 0;
        long regionRightWingRight = 0;
                
        ArrayList<SVGroup> variationList = new ArrayList();
        
        String indelType = "";
        if(inIndelType.equals("TD")){
            variationList = this.tandemList;
            regionFilterFileName=path.getParent().toString()+File.separator+sampleName+"_newReferenceTandem_regionFilter"+".txt";
        }else if(inIndelType.equals("ID")){
            variationList = this.indelList;
            regionFilterFileName=path.getParent().toString()+File.separator+sampleName+"_newReferenceIndel_regionFilter"+".txt";
        }else if(inIndelType.equals("IC")){
            variationList = this.intraTransList;
            regionFilterFileName=path.getParent().toString()+File.separator+sampleName+"_newReferenceIntraTrans_regionFilter"+".txt";
        }else if(inIndelType.equals("IT")){
            variationList = this.interTransList;
            regionFilterFileName=path.getParent().toString()+File.separator+sampleName+"_newReferencenterTrans_regionFilter"+".txt";
        }
        
        writer2 = new FileWriter(regionFilterFileName);
        
        /**
         * Start loop in array of sort coverage list
         */
        for(int i=0;i<variationList.size();i++){
            if(variationList.size() == 0){
                break;
            }
            
            SVGroup dummySVGroup = variationList.get(i);
            
            if(dummySVGroup.getNumCoverage() >= minPickCoverage){
                
                if(inIndelType.equals("TD")){
                    filename = path.getParent().toString()+File.separator+sampleName+"_junc"+i+"_newReferenceTandem.fa";                    
                }else if(inIndelType.equals("ID")){                   
                    filename = path.getParent().toString()+File.separator+sampleName+"_junc"+i+"_newReferenceIndel.fa";                    
                }else if(inIndelType.equals("IC")){                    
                    filename = path.getParent().toString()+File.separator+sampleName+"_junc"+i+"_newReferenceIntraTrans.fa";                    
                }else if(inIndelType.equals("IT")){                    
                    filename = path.getParent().toString()+File.separator+sampleName+"_junc"+i+"_newReferenceInterTrans.fa";                    
                }
                writer = new FileWriter(filename);
                /**
                 * pick first variation object to check for indel type and get 
                 */
                VariationV2 dummyVariationV2 = dummySVGroup.getVarList().get(0);        // pick first read of this group
                
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
                
                
//                long regionLeftWingLeft = startLeftWing-(extendSize*2);
//                if(startLeftWing<(extendSize*2)){
//                    regionLeftWingLeft = 0;
//                }
//                long regionLeftWingRight = startLeftWing+(extendSize*2);
//                
//                String region = chrNameF+":"+startLeftWing+
//                writer2.write(chrNameF);
                
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
                    
//                    regionLeftWingLeft = dummyVariationV2.getBreakpointF();
//                    regionLeftWingRight = startLeftWing-(extendSize*2);
//                    if(startLeftWing<(extendSize*2)){
//                        regionLeftWingRight = 0;
//                    }
//                    regionRightWingLeft = startRightWing+(extendSize*2);
//                    regionRightWingRight = dummyVariationV2.getBreakpointB();                    
                    
                }else if(strand.equals("01")){
                    String rightWingInvSeq = SequenceUtil.inverseSequence(rightWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String rightWingCompSeq = SequenceUtil.createComplimentV2(rightWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    newRefSeq = leftWingSeq + readSeq + rightWingCompSeq;
                    
//                    regionLeftWingLeft = startLeftWing-(extendSize*2);
//                    if(startLeftWing<(extendSize*2)){
//                        regionLeftWingLeft = 0;
//                    }
//                    regionLeftWingRight = dummyVariationV2.getBreakpointF();                    
//                    regionRightWingLeft = startRightWing+(extendSize*2);
//                    regionRightWingRight = dummyVariationV2.getBreakpointB();
                    
                }else if(strand.equals("10")){
                    String leftWingInvSeq = SequenceUtil.inverseSequence(leftWingSeq);                                // Do invert sequence (ATCG => GCTA)
                    String leftWingCompSeq = SequenceUtil.createComplimentV2(leftWingInvSeq);                       // Do compliment on invert sequence (GCTA => CGAT)

                    newRefSeq = leftWingCompSeq + readSeq + rightWingSeq;
                    
//                    regionLeftWingLeft = dummyVariationV2.getBreakpointF();
//                    regionLeftWingRight = startLeftWing-(extendSize*2);
//                    if(startLeftWing<(extendSize*2)){
//                        regionLeftWingRight = 0;
//                    }                    
//                    regionRightWingLeft = dummyVariationV2.getBreakpointB();
//                    regionRightWingRight = startRightWing+(extendSize*2);
                }

                String nameRef = "group"+groupCount+"_"+dummyVariationV2.readID;

                writer.write(">"+nameRef);
                writer.write("\n");
                writer.write(newRefSeq);
                writer.write("\n");
                
                /**
                 * prepare region data
                 * and write data
                 */
                regionLeftWingLeft = startLeftWing-(extendSize*2);
                if(startLeftWing<(extendSize*2)){
                    regionLeftWingLeft = 0;
                }
                regionLeftWingRight = dummyVariationV2.getBreakpointF();                    
                regionRightWingLeft = dummyVariationV2.getBreakpointB();
                regionRightWingRight = startRightWing+(extendSize*2);
                
                writer2.write(chrNameF+":"+regionLeftWingLeft+"-"+regionLeftWingRight+" "+chrNameB+":"+regionRightWingLeft+"-"+regionRightWingRight+"\n");
                /****************************************/
            }else{
                break;
            }
            writer.flush();
            writer.close();
        }
        writer2.flush();
        writer2.close();
        rbRef.close();
        
    }
    
    public ArrayList<String> createReferenceFromPreciseSVType(SVGroup inSVGroup, String refFile, String refIdxFile, int inExtendSize) throws IOException{
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
         * Implement save single junction file
         * 
         * ##Not finish
         */
        ArrayList<String> res = new ArrayList();
        RandomAccessFile rbRef = new RandomAccessFile(refFile,"r");
        
        RefFaIndex refFaIdx = new RefFaIndex(refIdxFile);

        int extendSize = inExtendSize;
        int groupCount = 0;


        /**
         * pick first variation object to check for indel type and get 
         */
        VariationV2 dummyVariationV2 = inSVGroup.getVarList().get(0);        // pick first read of this group

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

        long numLineRead = ((extendSize+unmatchFront)/refFaIdx.getLineBase(chrNameF))*3;
        if(numLineRead == 0) numLineRead=1*3;
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
        
        numLineRead = ((extendSize+unmatchBack)/refFaIdx.getLineBase(chrNameB))*3;
        if(numLineRead == 0) numLineRead=1*3;
        
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

        String nameRef = ""+dummyVariationV2.readID;
//        String newRefNumBaseBeforeBPF = "" + (leftWingSeq.length() + (dummyVariationV2.getBreakpointIndexF() + 1));
        String newRefNumBaseBeforeBPF = "" + ((leftWingSeq.length()-unmatchFront) + (dummyVariationV2.getBreakpointIndexF() + 1));
//        if(newRefNumBaseBeforeBPF.equals("263")){
//                            System.out.println();
//                        }
        res.add(nameRef);
        res.add(newRefSeq);
        res.add(newRefNumBaseBeforeBPF);
        
        rbRef.close();
        return res;
        
    }
    
    public void writeVisualizePreciseSVType(String readFile, String refFile, String refIdxFile, String SVType, int minPickCoverage, int inExtendSize) throws IOException{
        /**
         * this function will write the precise SVType report with single base resolution
         * user must define SV type Code which has 5 SV Type
         * TD = tandem, D = deletion, IA = intraTrans, IE = interTrans, CH = chimeric
         */
        
        Path path = Paths.get(readFile);
        String fileName = path.getFileName().toString();
        String[] dummy2 = fileName.split("_unmap");         // For more generic this should be split by . (and whole protocal should use . as a filed peaparation like A.unmap.emdup.bam So, we can get read name just split by .)
        String sampleName = dummy2[0];
        String filename = "";
        FileWriter writer;
        FileWriter writer2;
        int extendSize = inExtendSize;
        int groupCount = 0;
        ArrayList<SVGroup> variationList = new ArrayList();
        Map<ArrayList<SVGroup>,Boolean> variationMap = new LinkedHashMap();
        
        String indelType = "";
        if(SVType.equals("TD")){
            variationList = this.tandemList;
            filename = path.getParent().toString()+File.separator+sampleName+".Tandem.report";            
        }else if(SVType.equals("D")){
            variationList = this.deletionList;
            filename = path.getParent().toString()+File.separator+sampleName+".Deletion.report";            
        }else if(SVType.equals("IA")){
            variationMap = this.intraTransPairList;
            filename = path.getParent().toString()+File.separator+sampleName+".IntraTranslocation.report";            
        }else if(SVType.equals("IE")){
            variationMap = this.interTransPairList;
            filename = path.getParent().toString()+File.separator+sampleName+".InterTranslocation.report";            
        }else if(SVType.equals("CH")){
            variationList = this.chimericList;
            filename = path.getParent().toString()+File.separator+sampleName+".Chimeric.report";            
        }
        
        writer = new FileWriter(filename);
        /**
         * For tandem or deletion
         */
        if(SVType.equals("TD")||SVType.equals("D")){
            
            for(int i=0;i<variationList.size();i++){
//                if(i==12){
//                            System.out.println();
//                        }
                SVGroup dummySVGroup = variationList.get(i);
                ArrayList<String> res = createReferenceFromPreciseSVType(dummySVGroup,refFile,refIdxFile,extendSize);
                String refName = res.get(0);
                int numBaseBeforeNewRefBPF = Integer.parseInt(res.get(2));
                
                if(SVType.equals("TD")){
                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortTandemSummary()+"\n");
                }else if(SVType.equals("D")){
                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortDeletionSummary()+"\n");
                }

                ArrayList<VariationV2> varList = dummySVGroup.getVarList();
                
                StringBuilder newRef = new StringBuilder(res.get(1));
                int newRefBPIdxF = 0;
                int newRefBPIdxB = 0;
                
                for(int k=0;k<varList.size();k++){
                    VariationV2 var = varList.get(k);

                    String read = var.getReadSeq();
                    int numBaseBeforeReadBPF = var.getBreakpointIndexF()+1;
                    int numBaseBeforeReadBPB = var.getBreakpointIndexB();               // special for breakpoint back no need to plus 1 to get num base
                    int numAppend = numBaseBeforeNewRefBPF - numBaseBeforeReadBPF;
//                    int numAppend = extendSize;
                    
                    StringBuilder readBuilder = new StringBuilder(newRef.length());
                    for(int num=0;num<numAppend;num++){
                        readBuilder.append(" ");
                    }
                    readBuilder.append(read);
                    
                    if(k==0){
                        // calculate variable for reference
                        // prepare reference and write reference
                        newRefBPIdxF = numBaseBeforeNewRefBPF; //+ numBaseBeforeReadBPF;                        
                        newRef.insert(newRefBPIdxF, "|");
                        if(var.getBreakpointIndexB() - var.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB = newRefBPIdxF;
                            newRef.insert(newRefBPIdxF, "|");
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB = numAppend + numBaseBeforeReadBPB + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                            newRef.insert(newRefBPIdxB, "|");             
                        }
                        writer.write(newRef.toString());
                        writer.write("\n");
                    }
                    // calculate variable for read
                    newRefBPIdxF = numAppend + numBaseBeforeReadBPF;
                    if(var.getBreakpointIndexB() - var.getBreakpointIndexF() == 1){
                        // BPF and BPB is the same location we can insert at the same index
                        newRefBPIdxB = newRefBPIdxF;                            
                    }else{
                        // Has unmatch base between BPF and BPB
                        newRefBPIdxB = numAppend + numBaseBeforeReadBPB + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                    }
                    
//                    if(newRefBPIdxF == 212){
//                        System.out.println();
//                    }
                    readBuilder.insert(newRefBPIdxF, "|");
                    readBuilder.insert(newRefBPIdxB, "|");
                    
                    writer.write(readBuilder.toString());
                    writer.write("\n");    
                } 
            }
            
        }else if(SVType.equals("IA")||SVType.equals("IE")){
            /**
             * For inter and intra Translocation
             */
            String dotSeparator = "..........";
            int count = 1;
            for(Map.Entry<ArrayList<SVGroup>,Boolean> entry : variationMap.entrySet()){
                SVGroup frontSVGroup = entry.getKey().get(0);
                SVGroup backSVGroup = entry.getKey().get(1);
                
                //Prepare newRef Front SVGroup
                ArrayList<String> resF = createReferenceFromPreciseSVType(frontSVGroup,refFile,refIdxFile,inExtendSize);
                String refNameF = resF.get(0);
//                String newRefF = resF.get(1);
                
                int numBaseBeforeNewRefBPF_F = Integer.parseInt(resF.get(2));
                /********/
                
                //Prepare newRef Back SVGroup
                ArrayList<String> resB = createReferenceFromPreciseSVType(backSVGroup,refFile,refIdxFile,inExtendSize);
                String refNameB = resB.get(0);
//                String newRefB = resB.get(1);
                
                int numBaseBeforeNewRefBPF_B = Integer.parseInt(resB.get(2));
                /********/
                
                writer.write(">"+count+"\t"+refNameF+"_"+refNameB+"\tFront : "+frontSVGroup.shortSummary()+"\tBack : "+backSVGroup.shortSummary()+"\n");
                count++;
//                writer.write(newRefF+".........."+newRefB);
//                writer.write("\n");
                
                //Loop varlist pick big list as a main loop
                ArrayList<VariationV2> varListF = frontSVGroup.getVarList();
                ArrayList<VariationV2> varListB = backSVGroup.getVarList();
                
                int maxSize = 0;
                if(varListF.size()>=varListB.size()){
                    maxSize = varListF.size();
                }else{
                    maxSize = varListB.size();
                }
                StringBuilder newRefF = new StringBuilder(resF.get(1).length()+2);      // create stringbuilder for newrefF with size +2 of the refF size (plus 2 because we add "|" 2 time to denote the breakpoint)
                int newRefBPIdxF_F = 0;
                int newRefBPIdxB_F = 0;
                newRefF.append(resF.get(1));
                int referenceFrontLen = resF.get(1).length();
//                int numBaseBeforeReadBPF_F = 0;
                
                StringBuilder newRefB = new StringBuilder(resB.get(1).length()+2);       // create stringbuilder for newrefB with size +2 of the refB size (plus 2 because we add "|" 2 time to denote the breakpoint)
                int newRefBPIdxF_B = 0;
                int newRefBPIdxB_B = 0;
                newRefB.append(resB.get(1));
//                int numBaseBeforeReadBPF_B = 0;
                int remainBaseOnBackOfFront = 0;
                for(int k=0;k<maxSize;k++){
                    int numAppendF = 0;
                    int numAppendB = 0;
                    
                    StringBuilder readBuilder = new StringBuilder(newRefF.length()+10+newRefB.length());
                    StringBuilder readBuilderF = new StringBuilder(newRefF.capacity());
                    StringBuilder readBuilderB = new StringBuilder(newRefB.capacity());
                    
                    //StringBuilder must be separate to readF readB then combine at the end to make it more easy to to under stand and cut the sapeend part and complex calculate off
                    // Front
                    if(k < varListF.size()){
                        VariationV2 varF = varListF.get(k);
                        String readF = varF.getReadSeq();
                        int numBaseBeforeReadBPF_F = varF.getBreakpointIndexF()+1;
                        int numBaseBeforeReadBPB_F = varF.getBreakpointIndexB();
                        numAppendF = numBaseBeforeNewRefBPF_F - numBaseBeforeReadBPF_F;


                        for(int num=0;num<numAppendF;num++){
                            readBuilderF.append(" ");
                        }
                        readBuilderF.append(readF);
                        remainBaseOnBackOfFront = referenceFrontLen - (numAppendF + readF.length());
/********/                        
                        if(k==0){
                            // calculate variable for ref
                            // prepare reference and write reference
                            // focus on front
//                            if(numBaseBeforeNewRefBPF_F == 120){
//                                System.out.println();
//                            }
                            newRefBPIdxF_F = numBaseBeforeNewRefBPF_F; //+ numBaseBeforeReadBPF_F; 
                            newRefF.insert(newRefBPIdxF_F, "|");
                            if(varF.getBreakpointIndexB() - varF.getBreakpointIndexF() == 1){
                                // BPF and BPB is the same location we can insert at the same index
                                newRefBPIdxB_F = newRefBPIdxF_F;
                                newRefF.insert(newRefBPIdxF_F, "|");
                            }else{
                                // Has unmatch base between BPF and BPB
                                newRefBPIdxB_F = numAppendF + numBaseBeforeReadBPB_F + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                                newRefF.insert(newRefBPIdxB_F, "|");             
                            }
//                            writer.write(newRefF.toString());
//                            writer.write("\n");
                        }
                        
                        // calculate variable for read
                        newRefBPIdxF_F = numAppendF + numBaseBeforeReadBPF_F;
                        if(varF.getBreakpointIndexB() - varF.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB_F = newRefBPIdxF_F;                            
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB_F = numAppendF + numBaseBeforeReadBPB_F + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                        }
                        /*********/
                    }

                    //Back
                    if(k < varListB.size()){
                        VariationV2 varB = varListB.get(k);
                        String readB = varB.getReadSeq();
                        int numBaseBeforeReadBPF_B = varB.getBreakpointIndexF()+1; 
                        int numBaseBeforeReadBPB_B = varB.getBreakpointIndexB();
                        if(readBuilderF.toString().length()==0){
                            numAppendB = newRefF.length() + dotSeparator.length() + (numBaseBeforeNewRefBPF_B - numBaseBeforeReadBPF_B);     // plus 10 in this line came from number of dot that we add we add 10 dot to separate newRefF and newRefB
                        }else{
                            
                            numAppendB = remainBaseOnBackOfFront + dotSeparator.length() + (numBaseBeforeNewRefBPF_B - numBaseBeforeReadBPF_B);     // plus 10 in this line came from number of dot that we add we add 10 dot to separate newRefF and newRefB
                        }
                             
                        for(int num=0;num<numAppendB;num++){
                            readBuilderB.append(" ");
                        }
                        
                        readBuilderB.append(readB);
 /***********/                       
                        if(k==0){
                            // prepare reference and write reference
                            // focus on front
                            newRefBPIdxF_B = numBaseBeforeNewRefBPF_B; // + numBaseBeforeReadBPF_B; 
                            newRefB.insert(newRefBPIdxF_B, "|");
                            if(varB.getBreakpointIndexB() - varB.getBreakpointIndexF() == 1){
                                // BPF and BPB is the same location we can insert at the same index
                                newRefBPIdxB_B = newRefBPIdxF_B;
                                newRefB.insert(newRefBPIdxF_B, "|");
                            }else{
                                // Has unmatch base between BPF and BPB
                                newRefBPIdxB_B = (numAppendB - (remainBaseOnBackOfFront + dotSeparator.length())) + numBaseBeforeReadBPB_B + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                                newRefB.insert(newRefBPIdxB_B, "|");             
                            }
                            
                            writer.write(newRefF+dotSeparator+newRefB);
                            writer.write("\n");
                        }
                       
                        // calculate variable for read
                        newRefBPIdxF_B = numAppendB + numBaseBeforeReadBPF_B;
                        if(varB.getBreakpointIndexB() - varB.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB_B = newRefBPIdxF_B;                            
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB_B = numAppendB + numBaseBeforeReadBPB_B + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                        }
                        /************/
                    }
                    
                    // Fill "|" into readBuilderF and readBuilderB
                    if(readBuilderF.toString().length()!=0){
                        readBuilderF.insert(newRefBPIdxF_F, "|");
                        readBuilderF.insert(newRefBPIdxB_F, "|");
                    }else{
//                        int numToAppend = newRefF.toString().length();
//                        for(int r=0;r<numToAppend;r++){
//                            readBuilderF.append(" ");
//                        }   
                    }
                    
                    if(readBuilderB.toString().length()!=0){
                        readBuilderB.insert(newRefBPIdxF_B, "|");
                        readBuilderB.insert(newRefBPIdxB_B, "|");
                    }else{
//                        int numToAppend = newRefB.toString().length();
//                        for(int r=0;r<numToAppend;r++){
//                            readBuilderB.append(" ");
//                        }
                    }
//                    readBuilderF.insert(newRefBPIdxF_F, "|");
//                    readBuilderF.insert(newRefBPIdxB_F, "|");
//                    readBuilderB.insert(newRefBPIdxF_B, "|");
//                    readBuilderB.insert(newRefBPIdxB_B, "|");
                    
                    writer.write(readBuilderF.toString()+readBuilderB.toString());
                    writer.write("\n");
                }
            }
            
//            for(int i=0;i<variationList.size();i++){
//                SVGroup dummySVGroup = variationList.get(i);
//                ArrayList<String> res = createReferenceFromPreciseSVType(dummySVGroup,refFile,refIdxFile,inExtendSize);
//                String refName = res.get(0);
//                String newRef = res.get(1);
//                int numBaseBeforeNewRefBPF = Integer.parseInt(res.get(3));
//                
//                if(SVType.equals("TD")){
//                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortTandemSummary()+"\n");
//                }else if(SVType.equals("D")){
//                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortDeletionSummary()+"\n");
//                }
//                
//                writer.write(newRef);
//                writer.write("\n");
//                
//                ArrayList<VariationV2> varList = dummySVGroup.getVarList();
//                for(int k=0;k<varList.size();k++){
//                    VariationV2 var = varList.get(k);
//                    String read = var.getReadSeq();
//                    int numBaseBeforeReadBPF = var.getBreakpointIndexF()+1;
//                    int numAppend = numBaseBeforeNewRefBPF + numBaseBeforeReadBPF;
//                    
//                    StringBuilder readBuilder = new StringBuilder();                    
//                    for(int num=0;num<numAppend;num++){
//                        readBuilder.append(" ");
//                    }
//                    readBuilder.append(read);
//                    
//                    writer.write(readBuilder.toString());
//                    writer.write("\n");    
//                } 
//            }
        }else if(SVType.equals("CH")){
            
        }
        
        writer.flush();
        writer.close();
    }
    
    public void writeVisualizePreciseSVTypeWithAnnotation(String readFile, String refFile, String refFaIndex, String gffFile, String SVType, int minPickCoverage, int inExtendSize, int merSize) throws IOException{
        /**
         * this function will write the precise SVType report with single base resolution
         * user must define SV type Code which has 5 SV Type
         * TD = tandem, D = deletion, IA = intraIns, IE = interIns, CH = chimeric
         */
        
        /**
        * Read annotation file (GFF)
        */
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV4(gffFile,refFaIndex, "gene" , merSize);
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        /****************************/
        
        Path path = Paths.get(readFile);
        String fileName = path.getFileName().toString();
        String[] dummy2 = fileName.split("_unmap");         // For more generic this should be split by . (and whole protocal should use . as a filed peaparation like A.unmap.emdup.bam So, we can get read name just split by .)
        String sampleName = dummy2[0];
        String filename = "";
        String readmeFile = "";
        FileWriter writer;
        FileWriter readmeWriter;
        FileWriter writer2;
        int extendSize = inExtendSize;
        int groupCount = 0;
        ArrayList<SVGroup> variationList = new ArrayList();
        ArrayList<SVGroupPair> variationPairList = new ArrayList();
        
        String indelType = "";
        if(SVType.equals("TD")){
            variationList = this.tandemList;
            filename = path.getParent().toString()+File.separator+"Tandem"+File.separator+sampleName+".Tandem.Anno.report";
            readmeFile = path.getParent().toString()+File.separator+"Tandem"+File.separator+sampleName+".Tandem.Anno.readme.txt";
            File file = new File(filename);
            if(!file.getParentFile().exists()){
                try{
                    file.getParentFile().mkdirs();
                } 
                catch(SecurityException se){
                    System.out.println("Can not create directory : " + file.getParent());
                }
            }
        }else if(SVType.equals("D")){
            variationList = this.deletionList;
            filename = path.getParent().toString()+File.separator+"Deletion"+File.separator+sampleName+".Deletion.Anno.report";
            readmeFile = path.getParent().toString()+File.separator+"Deletion"+File.separator+sampleName+".Deletion.Anno.readme.txt";
            File file = new File(filename);
            if(!file.getParentFile().exists()){
                try{
                    file.getParentFile().mkdirs();
                } 
                catch(SecurityException se){
                    System.out.println("Can not create directory : " + file.getParent());
                }
            }
        }else if(SVType.equals("IA")){
            variationPairList = this.intraInsertionList;
            filename = path.getParent().toString()+File.separator+"IntraInsertion"+File.separator+sampleName+".IntraInsertion.Anno.report";
            readmeFile = path.getParent().toString()+File.separator+"IntraInsertion"+File.separator+sampleName+".IntraInsertion.Anno.readme.txt";
            File file = new File(filename);
            if(!file.getParentFile().exists()){
                try{
                    file.getParentFile().mkdirs();
                } 
                catch(SecurityException se){
                    System.out.println("Can not create directory : " + file.getParent());
                }
            }
        }else if(SVType.equals("IE")){
            variationPairList = this.interInsertionList;
            filename = path.getParent().toString()+File.separator+"InterInsertion"+File.separator+sampleName+".InterInsertion.Anno.report";
            readmeFile = path.getParent().toString()+File.separator+"InterInsertion"+File.separator+sampleName+".InterInsertion.Anno.readme.txt";
            File file = new File(filename);
            if(!file.getParentFile().exists()){
                try{
                    file.getParentFile().mkdirs();
                } 
                catch(SecurityException se){
                    System.out.println("Can not create directory : " + file.getParent());
                }
            }            
        }else if(SVType.equals("CH")){
            variationList = this.chimericList;
            filename = path.getParent().toString()+File.separator+"Chimeric"+File.separator+sampleName+".Chimeric.Anno.report";
            readmeFile = path.getParent().toString()+File.separator+"Chimeric"+File.separator+sampleName+".Chimeric.Anno.readme.txt";
            File file = new File(filename);
            if(!file.getParentFile().exists()){
                try{
                    file.getParentFile().mkdirs();
                } 
                catch(SecurityException se){
                    System.out.println("Can not create directory : " + file.getParent());
                }
            }
        }
        
        writer = new FileWriter(filename);
        readmeWriter = new FileWriter(readmeFile);
        readmeWriter.write(generateReadme(SVType));
        /**
         * For tandem or deletion
         */
        if(SVType.equals("TD")||SVType.equals("D")||SVType.equals("CH")){
            int count = 1;
            for(int i=0;i<variationList.size();i++){
//                if(i==12){
//                            System.out.println();
//                        }
                SVGroup dummySVGroup = variationList.get(i);
                ArrayList<String> res = createReferenceFromPreciseSVType(dummySVGroup,refFile,refFaIndex,extendSize);
                String refName = res.get(0);
                int numBaseBeforeNewRefBPF = Integer.parseInt(res.get(2));
                
                /**
                 * Annotation part
                 */
                ArrayList<String> annotation = getAnnotationOfSVGroupV3(dummySVGroup, refAnno, refAnnoIndex, "gene",merSize);
                String annoF = annotation.get(0);
                String annoB = annotation.get(1);
                /******************/
                
                if(SVType.equals("TD")){
                    writer.write(">"+count+"\t"+refName+"\t"+dummySVGroup.shortTandemSummary()+"\tAnnotation front : " + annoF + "\tAnnotation back : "+ annoB + "\n");                    
                }else if(SVType.equals("D")){
                    writer.write(">"+count+"\t"+refName+"\t"+dummySVGroup.shortDeletionSummary()+"\tAnnotation front : " + annoF + "\tAnnotation back : "+ annoB+"\n");
                }else if(SVType.equals("CH")){
                    writer.write(">"+count+"\t"+refName+"\t"+dummySVGroup.shortSummary()+"\tAnnotation front : " + annoF + "\tAnnotation back : "+ annoB+"\n");
                }
                count++;
                ArrayList<VariationV2> varList = dummySVGroup.getVarList();
                
                StringBuilder newRef = new StringBuilder(res.get(1));
                int newRefBPIdxF = 0;
                int newRefBPIdxB = 0;
                
                for(int k=0;k<varList.size();k++){
                    VariationV2 var = varList.get(k);

                    String read = var.getReadSeq();
                    int numBaseBeforeReadBPF = var.getBreakpointIndexF()+1;
                    int numBaseBeforeReadBPB = var.getBreakpointIndexB();               // special for breakpoint back no need to plus 1 to get num base
                    int numAppend = numBaseBeforeNewRefBPF - numBaseBeforeReadBPF;
//                    int numAppend = extendSize;
                    
                    StringBuilder readBuilder = new StringBuilder(newRef.length());
                    for(int num=0;num<numAppend;num++){
                        readBuilder.append(" ");
                    }
                    readBuilder.append(read);
                    
                    if(k==0){
                        // calculate variable for reference
                        // prepare reference and write reference
                        newRefBPIdxF = numBaseBeforeNewRefBPF; //+ numBaseBeforeReadBPF;                        
                        newRef.insert(newRefBPIdxF, "|");
                        if(var.getBreakpointIndexB() - var.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB = newRefBPIdxF;
                            newRef.insert(newRefBPIdxF, "|");
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB = numAppend + numBaseBeforeReadBPB + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                            newRef.insert(newRefBPIdxB, "|");             
                        }
                        writer.write(newRef.toString());
                        writer.write("\n");
                    }
                    // calculate variable for read
                    newRefBPIdxF = numAppend + numBaseBeforeReadBPF;
                    if(var.getBreakpointIndexB() - var.getBreakpointIndexF() == 1){
                        // BPF and BPB is the same location we can insert at the same index
                        newRefBPIdxB = newRefBPIdxF;                            
                    }else{
                        // Has unmatch base between BPF and BPB
                        newRefBPIdxB = numAppend + numBaseBeforeReadBPB + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                    }
                    
//                    if(newRefBPIdxF == 212){
//                        System.out.println();
//                    }
                    readBuilder.insert(newRefBPIdxF, "|");
                    readBuilder.insert(newRefBPIdxB, "|");
                    
                    writer.write(readBuilder.toString());
                    writer.write("\n");    
                }
                writer.write("\n");
            }
            
        }else if(SVType.equals("IA")||SVType.equals("IE")){
            /**
             * For inter and intra Translocation
             */
            String dotSeparator = "..........";
            int count = 1;
//            for(Map.Entry<ArrayList<SVGroup>,Boolean> entry : variationMap.entrySet()){
//                SVGroup frontSVGroup = entry.getKey().get(0);
//                SVGroup backSVGroup = entry.getKey().get(1);
                
            for(int i=0;i<variationPairList.size();i++){
                SVGroupPair svPair = variationPairList.get(i);
                SVGroup frontSVGroup = svPair.getFrontSVGroup();
                SVGroup backSVGroup = svPair.getBackSVGroup();
                
                //Prepare newRef Front SVGroup
                ArrayList<String> resF = createReferenceFromPreciseSVType(frontSVGroup,refFile,refFaIndex,inExtendSize);
                String refNameF = resF.get(0);
//                String newRefF = resF.get(1);
                
                int numBaseBeforeNewRefBPF_F = Integer.parseInt(resF.get(2));
                /********/
                
                //Prepare newRef Back SVGroup
                ArrayList<String> resB = createReferenceFromPreciseSVType(backSVGroup,refFile,refFaIndex,inExtendSize);
                String refNameB = resB.get(0);
//                String newRefB = resB.get(1);
                
                int numBaseBeforeNewRefBPF_B = Integer.parseInt(resB.get(2));
                /********/
                
                /**
                 * Annotation part
                 */
                ArrayList<String> annotationF = getAnnotationOfSVGroupV3(frontSVGroup, refAnno, refAnnoIndex, "gene", merSize);
                String annoF_F = annotationF.get(0);
                String annoB_F = annotationF.get(1);
                ArrayList<String> annotationB = getAnnotationOfSVGroupV3(backSVGroup, refAnno, refAnnoIndex, "gene", merSize);
                String annoF_B = annotationB.get(0);
                String annoB_B = annotationB.get(1);
                /******************/
                
                
                writer.write(">"+count+"\t"+refNameF+"_"+refNameB+"\tFront : "+frontSVGroup.shortSummary()+"\tAnnotation Front : "+annoF_F+"\tAnnotation back : "+annoB_F+"\t||\tBack : "+backSVGroup.shortSummary()+"\tAnnotation Front : "+annoF_B+"\tAnnotation back : "+annoB_B
                        +"\t"+svPair.getInsertionSize()+"\t"+svPair.getInsertionJunction()+"\n");
                count++;
//                writer.write(newRefF+".........."+newRefB);
//                writer.write("\n");
                
                //Loop varlist pick big list as a main loop
                ArrayList<VariationV2> varListF = frontSVGroup.getVarList();
                ArrayList<VariationV2> varListB = backSVGroup.getVarList();
                
                int maxSize = 0;
                if(varListF.size()>=varListB.size()){
                    maxSize = varListF.size();
                }else{
                    maxSize = varListB.size();
                }
                StringBuilder newRefF = new StringBuilder(resF.get(1).length()+2);      // create stringbuilder for newrefF with size +2 of the refF size (plus 2 because we add "|" 2 time to denote the breakpoint)
                int newRefBPIdxF_F = 0;
                int newRefBPIdxB_F = 0;
                newRefF.append(resF.get(1));
                int referenceFrontLen = resF.get(1).length();
//                int numBaseBeforeReadBPF_F = 0;
                
                StringBuilder newRefB = new StringBuilder(resB.get(1).length()+2);       // create stringbuilder for newrefB with size +2 of the refB size (plus 2 because we add "|" 2 time to denote the breakpoint)
                int newRefBPIdxF_B = 0;
                int newRefBPIdxB_B = 0;
                newRefB.append(resB.get(1));
//                int numBaseBeforeReadBPF_B = 0;
                int remainBaseOnBackOfFront = 0;
                for(int k=0;k<maxSize;k++){
                    int numAppendF = 0;
                    int numAppendB = 0;
                    
                    StringBuilder readBuilder = new StringBuilder(newRefF.length()+10+newRefB.length());
                    StringBuilder readBuilderF = new StringBuilder(newRefF.capacity());
                    StringBuilder readBuilderB = new StringBuilder(newRefB.capacity());
                    
                    //StringBuilder must be separate to readF readB then combine at the end to make it more easy to to under stand and cut the sapeend part and complex calculate off
                    // Front
                    if(k < varListF.size()){
                        VariationV2 varF = varListF.get(k);
                        String readF = varF.getReadSeq();
                        int numBaseBeforeReadBPF_F = varF.getBreakpointIndexF()+1;
                        int numBaseBeforeReadBPB_F = varF.getBreakpointIndexB();
                        numAppendF = numBaseBeforeNewRefBPF_F - numBaseBeforeReadBPF_F;


                        for(int num=0;num<numAppendF;num++){
                            readBuilderF.append(" ");
                        }
                        readBuilderF.append(readF);
                        remainBaseOnBackOfFront = referenceFrontLen - (numAppendF + readF.length());
/********/                        
                        if(k==0){
                            // calculate variable for ref
                            // prepare reference and write reference
                            // focus on front
//                            if(numBaseBeforeNewRefBPF_F == 120){
//                                System.out.println();
//                            }
                            newRefBPIdxF_F = numBaseBeforeNewRefBPF_F; //+ numBaseBeforeReadBPF_F; 
                            newRefF.insert(newRefBPIdxF_F, "|");
                            if(varF.getBreakpointIndexB() - varF.getBreakpointIndexF() == 1){
                                // BPF and BPB is the same location we can insert at the same index
                                newRefBPIdxB_F = newRefBPIdxF_F;
                                newRefF.insert(newRefBPIdxF_F, "|");
                            }else{
                                // Has unmatch base between BPF and BPB
                                newRefBPIdxB_F = numAppendF + numBaseBeforeReadBPB_F + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                                newRefF.insert(newRefBPIdxB_F, "|");             
                            }
//                            writer.write(newRefF.toString());
//                            writer.write("\n");
                        }
                        
                        // calculate variable for read
                        newRefBPIdxF_F = numAppendF + numBaseBeforeReadBPF_F;
                        if(varF.getBreakpointIndexB() - varF.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB_F = newRefBPIdxF_F;                            
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB_F = numAppendF + numBaseBeforeReadBPB_F + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                        }
                        /*********/
                    }

                    //Back
                    if(k < varListB.size()){
                        VariationV2 varB = varListB.get(k);
                        String readB = varB.getReadSeq();
                        int numBaseBeforeReadBPF_B = varB.getBreakpointIndexF()+1; 
                        int numBaseBeforeReadBPB_B = varB.getBreakpointIndexB();
                        if(readBuilderF.toString().length()==0){
                            numAppendB = newRefF.length() + dotSeparator.length() + (numBaseBeforeNewRefBPF_B - numBaseBeforeReadBPF_B);     // plus 10 in this line came from number of dot that we add we add 10 dot to separate newRefF and newRefB
                        }else{
                            
                            numAppendB = remainBaseOnBackOfFront + dotSeparator.length() + (numBaseBeforeNewRefBPF_B - numBaseBeforeReadBPF_B);     // plus 10 in this line came from number of dot that we add we add 10 dot to separate newRefF and newRefB
                        }
                             
                        for(int num=0;num<numAppendB;num++){
                            readBuilderB.append(" ");
                        }
                        
                        readBuilderB.append(readB);
 /***********/                       
                        if(k==0){
                            // prepare reference and write reference
                            // focus on front
                            newRefBPIdxF_B = numBaseBeforeNewRefBPF_B; // + numBaseBeforeReadBPF_B; 
                            newRefB.insert(newRefBPIdxF_B, "|");
                            if(varB.getBreakpointIndexB() - varB.getBreakpointIndexF() == 1){
                                // BPF and BPB is the same location we can insert at the same index
                                newRefBPIdxB_B = newRefBPIdxF_B;
                                newRefB.insert(newRefBPIdxF_B, "|");
                            }else{
                                // Has unmatch base between BPF and BPB
                                newRefBPIdxB_B = (numAppendB - (remainBaseOnBackOfFront + dotSeparator.length())) + numBaseBeforeReadBPB_B + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF
                                newRefB.insert(newRefBPIdxB_B, "|");             
                            }
                            
                            writer.write(newRefF+dotSeparator+newRefB);
                            writer.write("\n");
                        }
                       
                        // calculate variable for read
                        newRefBPIdxF_B = numAppendB + numBaseBeforeReadBPF_B;
                        if(varB.getBreakpointIndexB() - varB.getBreakpointIndexF() == 1){
                            // BPF and BPB is the same location we can insert at the same index
                            newRefBPIdxB_B = newRefBPIdxF_B;                            
                        }else{
                            // Has unmatch base between BPF and BPB
                            newRefBPIdxB_B = numAppendB + numBaseBeforeReadBPB_B + 1; // plus 1 because breakpint is far from each other there is some shift base when we add "|" on BPF                                       
                        }
                        /************/
                    }
                    
                    // Fill "|" into readBuilderF and readBuilderB
                    if(readBuilderF.toString().length()!=0){
                        readBuilderF.insert(newRefBPIdxF_F, "|");
                        readBuilderF.insert(newRefBPIdxB_F, "|");
                    }else{
//                        int numToAppend = newRefF.toString().length();
//                        for(int r=0;r<numToAppend;r++){
//                            readBuilderF.append(" ");
//                        }   
                    }
                    
                    if(readBuilderB.toString().length()!=0){
                        readBuilderB.insert(newRefBPIdxF_B, "|");
                        readBuilderB.insert(newRefBPIdxB_B, "|");
                    }else{
//                        int numToAppend = newRefB.toString().length();
//                        for(int r=0;r<numToAppend;r++){
//                            readBuilderB.append(" ");
//                        }
                    }
//                    readBuilderF.insert(newRefBPIdxF_F, "|");
//                    readBuilderF.insert(newRefBPIdxB_F, "|");
//                    readBuilderB.insert(newRefBPIdxF_B, "|");
//                    readBuilderB.insert(newRefBPIdxB_B, "|");
                    
                    writer.write(readBuilderF.toString()+readBuilderB.toString());
                    writer.write("\n");
                }
                writer.write("\n");
            }
            
//            for(int i=0;i<variationList.size();i++){
//                SVGroup dummySVGroup = variationList.get(i);
//                ArrayList<String> res = createReferenceFromPreciseSVType(dummySVGroup,refFile,refIdxFile,inExtendSize);
//                String refName = res.get(0);
//                String newRef = res.get(1);
//                int numBaseBeforeNewRefBPF = Integer.parseInt(res.get(3));
//                
//                if(SVType.equals("TD")){
//                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortTandemSummary()+"\n");
//                }else if(SVType.equals("D")){
//                    writer.write(">"+i+"\t"+refName+"\t"+dummySVGroup.shortDeletionSummary()+"\n");
//                }
//                
//                writer.write(newRef);
//                writer.write("\n");
//                
//                ArrayList<VariationV2> varList = dummySVGroup.getVarList();
//                for(int k=0;k<varList.size();k++){
//                    VariationV2 var = varList.get(k);
//                    String read = var.getReadSeq();
//                    int numBaseBeforeReadBPF = var.getBreakpointIndexF()+1;
//                    int numAppend = numBaseBeforeNewRefBPF + numBaseBeforeReadBPF;
//                    
//                    StringBuilder readBuilder = new StringBuilder();                    
//                    for(int num=0;num<numAppend;num++){
//                        readBuilder.append(" ");
//                    }
//                    readBuilder.append(read);
//                    
//                    writer.write(readBuilder.toString());
//                    writer.write("\n");    
//                } 
//            }
        }
        
        writer.flush();
        writer.close();
        readmeWriter.flush();
        readmeWriter.close();
    }
    
    public ArrayList<String> getAnnotationOfSVGroup(SVGroup svGroup, ReferenceAnnotation refAnno, Map<Integer,Annotation> refAnnoIndex, String annotationField, int merSize) throws IOException{
            
//        /**
//        * Read annotation file (GFF)
//        */
//        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaIndex, annotationField , merSize);
//        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
//        /****************************/
        
        /**
         * annotation part
         * Find annotation for this svGroup
         */
        ArrayList<String> result = new ArrayList();
        if(!svGroup.isIdentityFlag()){
            svGroup.defineIdentity();
        }

        long avgIniPosF = svGroup.getRPF();
        long avgLastPosB = svGroup.getRPB();
        
        long strandChrF = ((long)svGroup.getStrandF()<<11)+(long)svGroup.getNumChrF();
        long strandChrB = ((long)svGroup.getStrandB()<<11)+(long)svGroup.getNumChrB();

        long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
        long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
        long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
        long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            
//        long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
//        long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//        long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
//        long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

        int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
        int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

        String strAnnoF = "null";
        String strAnnoB = "null";
        if(annoGroupIndexF >= 0){
            Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
            svGroup.setAnnoF(annoF);
            strAnnoF = annoF.toString();
        }
        if(annoGroupIndexB >= 0){
            Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
            svGroup.setAnnoB(annoB);
            strAnnoB = annoB.toString();
        }

        /**********************/
        result.add(strAnnoF);   // index 0 is front annotation
        result.add(strAnnoB);   // index 1 is back annotation

        return result;
    }
    
    public ArrayList<String> getAnnotationOfSVGroupV2(SVGroup svGroup, ReferenceAnnotation refAnno, Map<Integer,Annotation> refAnnoIndex, String annotationField, int merSize) throws IOException{
            
//        /**
//        * Read annotation file (GFF)
//        */
//        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaIndex, annotationField , merSize);
//        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
//        /****************************/
        
        /**
         * annotation part
         * Find annotation for this svGroup
         */
        boolean firstFoundFlagF = true;
        boolean firstFoundFlagB = true;
        ArrayList<String> strandPattern = new ArrayList();
        
        if(svGroup.isPpFlag()){
            strandPattern.add("++");
        }
        if(svGroup.isMmFlag()){
            strandPattern.add("--");
        }
        if(svGroup.isPmFlag()){
            strandPattern.add("+-");
        }
        if(svGroup.isMpFlag()){
            strandPattern.add("-+");
        }
        
        
        ArrayList<String> result = new ArrayList();
        if(!svGroup.isIdentityFlag()){
            svGroup.defineIdentity();
        }

        long strandF = 0;
        long strandB = 0;
        String strAnnoF = "null";
        String strAnnoB = "null";
        
        for(int i=0;i<strandPattern.size();i++){
            if(strandPattern.get(i).equals("++")){
                strandF = 0;
                strandB = 0;
            }else if(strandPattern.get(i).equals("--")){
                strandF = 1;
                strandB = 1;
            }else if(strandPattern.get(i).equals("+-")){
                strandF = 0;
                strandB = 1;
            }else if(strandPattern.get(i).equals("-+")){
                strandF = 1;
                strandB = 0;
            }
        
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long strandChrF = (strandF<<11)+(long)svGroup.getNumChrF();
            long strandChrB = (strandB<<11)+(long)svGroup.getNumChrB();

            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
            long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

    //        long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
    //        long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
    //        long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
    //        long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                if(firstFoundFlagF==true){
                    strAnnoF = annoF.toString();
                    firstFoundFlagF = false;
                }else{
                    strAnnoF = strAnnoF + " | " + annoF.toString();
                }               
            }
            
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                if(firstFoundFlagB==true){
                    strAnnoB = annoB.toString();
                    firstFoundFlagB = false;
                }else{
                    strAnnoB = strAnnoB + " | " + annoB.toString();
                }
            }

            /**********************/
            
        }
        result.add(strAnnoF);   // index 0 is front annotation
        result.add(strAnnoB);   // index 1 is back annotation
        return result;
    }
    
    public ArrayList<String> getAnnotationOfSVGroupV3(SVGroup svGroup, ReferenceAnnotation refAnno, Map<Integer,Annotation> refAnnoIndex, String annotationField, int merSize) throws IOException{
        /**
         * Annotated both strand + and - for every SVGroup even SVgroup it self did not have - strand or + strand
         */
//        /**
//        * Read annotation file (GFF) 
//        */
//        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaIndex, annotationField , merSize);
//        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
//        /****************************/
        
        /**
         * annotation part
         * Find annotation for this svGroup
         */
        boolean firstFoundFlagF = true;
        boolean firstFoundFlagB = true;
//        ArrayList<String> strandPattern = new ArrayList();
//        
//        if(svGroup.isPpFlag()){
//            strandPattern.add("++");
//        }
//        if(svGroup.isMmFlag()){
//            strandPattern.add("--");
//        }
//        if(svGroup.isPmFlag()){
//            strandPattern.add("+-");
//        }
//        if(svGroup.isMpFlag()){
//            strandPattern.add("-+");
//        }
        
        
        ArrayList<String> result = new ArrayList();
        if(!svGroup.isIdentityFlag()){
            svGroup.defineIdentity();
        }

        long strandF = 0;
        long strandB = 0;
        String strAnnoF = "null";
        String strAnnoB = "null";
        
        // loop 2 time fist stran + and second strand -
        for(int i=0;i<2;i++){
            if(i==0){
                strandF = 0;
                strandB = 0;
            }else if(i==1){
                strandF = 1;
                strandF = 1;
            }
            
//            if(strandPattern.get(i).equals("++")){
//                strandF = 0;
//                strandB = 0;
//            }else if(strandPattern.get(i).equals("--")){
//                strandF = 1;
//                strandB = 1;
//            }else if(strandPattern.get(i).equals("+-")){
//                strandF = 0;
//                strandB = 1;
//            }else if(strandPattern.get(i).equals("-+")){
//                strandF = 1;
//                strandB = 0;
//            }
        
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long strandChrF = (strandF<<11)+(long)svGroup.getNumChrF();
            long strandChrB = (strandB<<11)+(long)svGroup.getNumChrB();

            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
            long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

    //        long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
    //        long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
    //        long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
    //        long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                if(firstFoundFlagF==true){
                    strAnnoF = annoF.toString();
                    firstFoundFlagF = false;
                }else{
                    strAnnoF = strAnnoF + " | " + annoF.toString();
                }               
            }
            
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                if(firstFoundFlagB==true){
                    strAnnoB = annoB.toString();
                    firstFoundFlagB = false;
                }else{
                    strAnnoB = strAnnoB + " | " + annoB.toString();
                }
            }

            /**********************/
            
        }
        result.add(strAnnoF);   // index 0 is front annotation
        result.add(strAnnoB);   // index 1 is back annotation
        return result;
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
        Collections.sort(this.tandemList,SVGroup.CoverageComparator);
        Collections.sort(this.indelList,SVGroup.CoverageComparator);
        Collections.sort(this.intraTransList,SVGroup.CoverageComparator);
        Collections.sort(this.interTransList,SVGroup.CoverageComparator);        
    }
    
    public void classifyPreciseSVType(int coverageThreshold){
        /**
         * classify sv type of each sv candidate group
         * We hasve 5 precise type of SV (but this function will classify those SV type in more precisely)
         *  1. tandem
         *  2. indel
         *  3. intraIns
         *  4. interIns
         *  5. chimeric
         * 
         * Implement relation analysis algorithm for classify precise SV type
         * 
         * 
         */
        this.sameChrSVGroup = new ArrayList();
        this.diffChrSVGroup = new ArrayList();
//        int count=0;
        for(Map.Entry<String, Map<String,SVGroup>> entry1 : this.coverageMap.entrySet()){
            Map<String,SVGroup> coverageMapII = entry1.getValue();
            
            for(Map.Entry<String,SVGroup> entry2 : coverageMapII.entrySet()){
                SVGroup svGroup = entry2.getValue();
//                System.out.println(count++);
//                if(count==26){
//                    System.out.println();
//                }
                if(svGroup.getNumCoverage() >= coverageThreshold){
                    svGroup.defineIdentity();

                    if(svGroup.isDiffChrFlag()){
                        diffChrSVGroup.add(svGroup);
                    }else if(svGroup.isSameChrFlag()){
                        sameChrSVGroup.add(svGroup);
                    }
                }
                
                
                
//                if(svGroup.getSVType().equals("tandem")){
//                    this.tandemList.add(svGroup);
//                }else if(svGroup.getSVType().equals("indel")){
//                    this.indelList.add(svGroup);
//                }else if(svGroup.getSVType().equals("intraTrans")){
//                    this.intraTransList.add(svGroup);
//                }else if(svGroup.getSVType().equals("interTrans")){
//                    this.interTransList.add(svGroup);
//                }
            }
        }
        
        Collections.sort(diffChrSVGroup,SVGroup.CoverageComparator);
        Collections.sort(sameChrSVGroup,SVGroup.CoverageComparator);
        
        /**
         * Classifying precise SV type
         *  1. Start with calculate distance 
         *  2. classify same chromosome SV type
         *  3. classify different chromosome SV type
         */
        double eudistance = 0;
        for(int i=0;i<sameChrSVGroup.size();i++){
            SVGroup mainSV = sameChrSVGroup.get(i);
            for(int j=i;j<sameChrSVGroup.size();j++){
                SVGroup subSV = sameChrSVGroup.get(j);
                //calcualte 2D euclidean distance
//                double frontDiffSqrt = Math.pow((subSV.getFrontCode() - mainSV.getFrontCode()),2);
//                double backDiffSqrt = Math.pow((subSV.getBackCode() - mainSV.getBackCode()),2);
//                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
                /*******/
                if(mainSV.getEuCalFrontCode() == subSV.getEuCalFrontCode() && mainSV.getEuCalBackCode() == mainSV.getEuCalBackCode()){
                    eudistance = 0;
                }else{
                    //calcualte 2D reverse euclidean distance (switch subtraction order to x1-y2 and x2-y1)
                    double frontDiffSqrt = Math.pow((subSV.getEuCalBackCode() - mainSV.getEuCalFrontCode()),2);
                    double backDiffSqrt = Math.pow((subSV.getEuCalFrontCode() - mainSV.getEuCalBackCode()),2);
                    eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
                    /*******/
                }
                
                if(i==j){
                    //add eudistance one time if it is the same SVGroup
                    mainSV.addEuclideanDistance(eudistance);
                }else{
                    mainSV.addEuclideanDistance(eudistance);
                    subSV.addEuclideanDistance(eudistance);
                }
                 
            }    
        }
        
        for(int i=0;i<diffChrSVGroup.size();i++){
            SVGroup mainSV = diffChrSVGroup.get(i);
            for(int j=i;j<diffChrSVGroup.size();j++){
                SVGroup subSV = diffChrSVGroup.get(j);
                //calcualte 2D euclidean distance
//                double frontDiffSqrt = Math.pow((subSV.getFrontCode() - mainSV.getFrontCode()),2);
//                double backDiffSqrt = Math.pow((subSV.getBackCode() - mainSV.getBackCode()),2);
//                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
                /*******/
                
                //calcualte 2D reverse euclidean distance (switch subtraction order to x1-y2 and x2-y1)
                double frontDiffSqrt = Math.pow((subSV.getEuCalBackCode() - mainSV.getEuCalFrontCode()),2);
                double backDiffSqrt = Math.pow((subSV.getEuCalFrontCode() - mainSV.getEuCalBackCode()),2);
                eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
                /*******/
                
                if(i==j){
                    //add eudistance one time if it is the same SVGroup
                    mainSV.addEuclideanDistance(eudistance);
                }else{
                    mainSV.addEuclideanDistance(eudistance);
                    subSV.addEuclideanDistance(eudistance);
                }
                 
            }    
        }
        /**
         * identify same chromosome svtype
         */
        Map<Integer,Boolean> skipIndex = new LinkedHashMap();
        for(int i=0;i<this.sameChrSVGroup.size();i++){
//            if(skipIndex.containsKey(i)){
//                continue;
//            }
            SVGroup svGroupMain = this.sameChrSVGroup.get(i);
            ArrayList<Double> euList = svGroupMain.getSameChrEnclideanDistance();
            ArrayList<Double> dummyEuList = new ArrayList<Double>(euList);

            for(int k=0;k<this.sameChrSVGroup.size();k++){
                int minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                double minEuDis = dummyEuList.get(minIndex);
                
                if(minEuDis == 0){
                    // erase minindex if minEudis is 0 then repick
                    dummyEuList.remove(minIndex);
                    minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                    minEuDis = dummyEuList.get(minIndex);
                }
                
                int realMinIndex = euList.indexOf(minEuDis);    // true min index for write report
                SVGroup svGroupSub = this.sameChrSVGroup.get(realMinIndex);

                /**
                 * classify SV type
                 */

                String svType = identifySameChrPreciseSVType(svGroupMain,svGroupSub);  // this function will idetify SV type and put the SVGroup in to the correct SVtype list
                
                if(svType != null){
                    
                    
                    break;
                }
                
//                if(svType.equals("intraTrans")||svType.equals("interTrans")){
//                    /**
//                     * In case of inter and intra it will add SGroup into a list as a pair of SVGroup which compose of front and back svGroup. 
//                     * That's mean both mainSVGroup and sunSVGroup has been added
//                     * So, there is no reason to consider subSVGroup again when loop to it so we add the index of subSVGroup to skipIndex map
//                     * Checking the map and skip SVGroup that already consider when it come
//                     */
//                    skipIndex.put(realMinIndex, true);
//                }
                dummyEuList.remove(minIndex);
            }
        }
        
        /**
         * identify diff chromosome svtype
         */
        for(int i=0;i<this.diffChrSVGroup.size();i++){
//            if(skipIndex.containsKey(i)){
//                continue;
//            }
            SVGroup svGroupMain = this.diffChrSVGroup.get(i);
            ArrayList<Double> euList = svGroupMain.getDiffChrEuclideanDistance();
            ArrayList<Double> dummyEuList = new ArrayList<Double>(euList);

            for(int k=0;k<this.diffChrSVGroup.size();k++){
                int minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                double minEuDis = dummyEuList.get(minIndex);
                
                if(minEuDis == 0){
                    // erase minindex if minEudis is 0 then repick
                    dummyEuList.remove(minIndex);
                    minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                    minEuDis = dummyEuList.get(minIndex);
                }
                
                int realMinIndex = euList.indexOf(minEuDis);    // true min index for write report
                SVGroup svGroupSub = this.diffChrSVGroup.get(realMinIndex);

                /**
                 * classify SV type
                 */

                String svType = identifyDiffChrPreciseSVType(svGroupMain,svGroupSub);  // this function will idetify SV type and put the SVGroup in to the correct SVtype list
                
                if(svType != null){
                    
                    
                    break;
                }
                
//                if(svType.equals("intraTrans")||svType.equals("interTrans")){
//                    /**
//                     * In case of inter and intra it will add SGroup into a list as a pair of SVGroup which compose of front and back svGroup. 
//                     * That's mean both mainSVGroup and sunSVGroup has been added
//                     * So, there is no reason to consider subSVGroup again when loop to it so we add the index of subSVGroup to skipIndex map
//                     * Checking the map and skip SVGroup that already consider when it come
//                     */
//                    skipIndex.put(realMinIndex, true);
//                }
                dummyEuList.remove(minIndex);
            }
        }
        
        /*******************/
        
        
//        /**
//         * Analyze correctness of SVType (tandem and delete) for each svGroup 
//         * 
//         * Loop each svType list and check with our ideal hypothesis we call "relation of coverage with correctness"
//         * 
//         * Utilize samtools view command to get coverage of each specific region (our breakpoint)
//         */
//        
//        String baseCommand = "samtools view " + bamFile;
//        
//        // Loop deletion list
//        for(int i=0;i<this.deletionList.size();i++){
//            String addCommand = "";
//            SVGroup svGroup = deletionList.get(i);
//            
//            // query Coverage of frontbreakpoint
//            long breakpointF = svGroup.getRPF();
//            String chrNameF = svGroup.getChrF();
//            addCommand = "chr"+chrNameF+":"+breakpointF+"-"+breakpointF;
//            String command = baseCommand + " " + addCommand + " | wc";
//            
//            String resultF = executeCommandLine(command);
//            String[] splitResF = resultF.split("\\s+");
//            int coverageF = Integer.parseInt(splitResF[1]);
//            
//            // query Coverage f backreakpoint
//            long breakpointB = svGroup.getRPB();
//            String chrNameB = svGroup.getChrF();
//            addCommand = "chr"+chrNameB+":"+breakpointB+"-"+breakpointB;
//            command = baseCommand + " " + addCommand + " | wc";
//            
//            String resultB = executeCommandLine(command);
//            String[] splitResB = resultB.split("\\s+");
//            int coverageB = Integer.parseInt(splitResB[1]);
//            
//            // Add sum of two coverage to delectioncorrectness field if SVGroup
//            svGroup.setDeletionCorrectness(coverageF + coverageB);   
//        }
//        
//        
//        // Loop tandem List
//        for(int i=0;i<this.tandemList.size();i++){
//            String addCommand = "";
//            SVGroup svGroup = tandemList.get(i);
//            
//            // query Coverage of frontbreakpoint
//            long breakpointF = svGroup.getRPF();
//            String chrNameF = svGroup.getChrF();
//            addCommand = "chr"+chrNameF+":"+breakpointF+"-"+breakpointF;
//            String command = baseCommand + " " + addCommand + " | wc";
//            
//            String resultF = executeCommandLine(command);
//            String[] splitResF = resultF.split("\\s+");
//            int coverageF = Integer.parseInt(splitResF[1]);
//            
//            // query Coverage f backreakpoint
//            long breakpointB = svGroup.getRPB();
//            String chrNameB = svGroup.getChrF();
//            addCommand = "chr"+chrNameB+":"+breakpointB+"-"+breakpointB;
//            command = baseCommand + " " + addCommand + " | wc";
//            
//            String resultB = executeCommandLine(command);
//            String[] splitResB = resultB.split("\\s+");
//            int coverageB = Integer.parseInt(splitResB[1]);
//            
//            // Add sum of two coverage to tandemcorrectness field if SVGroup
//            svGroup.setTandemCorrectness(coverageF + coverageB);   
//        }
//        
//        Collections.sort(this.deletionList,SVGroup.CoverageCorrectnessDeletionComparator);
//        Collections.sort(this.tandemList,SVGroup.CoverageCorrectnessTandemComparator);
//        /********************/


        /**
         * Grouping reverse and normal in to the same Group 
         * calculate distant with code that compose of chromosome and breakpoint
         */
//        for(int i=0;i<sameChrSVGroup.size();i++){
//            SVGroup mainSV = sameChrSVGroup.get(i);
//            for(int j=i;j<sameChrSVGroup.size();j++){
//                SVGroup subSV = sameChrSVGroup.get(j);
//                //calcualte 2D euclidean distance
//                double frontDiffSqrt = Math.pow((subSV.getFrontCode() - mainSV.getFrontCode()),2);
//                double backDiffSqrt = Math.pow((subSV.getBackCode() - mainSV.getBackCode()),2);
//                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
//                /*******/
//                
//                //calcualte 2D reverse euclidean distance (switch subtraction order to x1-y2 and x2-y1)
////                double frontDiffSqrt = Math.pow((subSV.getEuCalBackCode() - mainSV.getEuCalFrontCode()),2);
////                double backDiffSqrt = Math.pow((subSV.getEuCalFrontCode() - mainSV.getEuCalBackCode()),2);
////                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
//                /*******/
//                
//                if(i==j){
//                    //add eudistance one time if it is the same SVGroup
//                    mainSV.addEuclideanDistance(eudistance);
//                }else{
//                    mainSV.addEuclideanDistance(eudistance);
//                    subSV.addEuclideanDistance(eudistance);
//                }
//                 
//            }    
//        }
//        
//        for(int i=0;i<diffChrSVGroup.size();i++){
//            SVGroup mainSV = diffChrSVGroup.get(i);
//            for(int j=i;j<diffChrSVGroup.size();j++){
//                SVGroup subSV = diffChrSVGroup.get(j);
//                //calcualte 2D euclidean distance
//                double frontDiffSqrt = Math.pow((subSV.getFrontCode() - mainSV.getFrontCode()),2);
//                double backDiffSqrt = Math.pow((subSV.getBackCode() - mainSV.getBackCode()),2);
//                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
//                /*******/
//                
//                //calcualte 2D reverse euclidean distance (switch subtraction order to x1-y2 and x2-y1)
////                double frontDiffSqrt = Math.pow((subSV.getEuCalBackCode() - mainSV.getEuCalFrontCode()),2);
////                double backDiffSqrt = Math.pow((subSV.getEuCalFrontCode() - mainSV.getEuCalBackCode()),2);
////                double eudistance = Math.sqrt((frontDiffSqrt + backDiffSqrt));
//                /*******/
//                
//                if(i==j){
//                    //add eudistance one time if it is the same SVGroup
//                    mainSV.addEuclideanDistance(eudistance);
//                }else{
//                    mainSV.addEuclideanDistance(eudistance);
//                    subSV.addEuclideanDistance(eudistance);
//                }
//                 
//            }    
//        }
        
        
        /********************/
        
        
        // sort list of SVGroup by number of coverage fron high to low (descending) comparable has benn inplement in SVGroup object
//        Collections.sort(this.tandemList,SVGroup.CoverageComparator);
//        Collections.sort(this.indelList,SVGroup.CoverageComparator);
//        Collections.sort(this.intraTransList,SVGroup.CoverageComparator);
//        Collections.sort(this.interTransList,SVGroup.CoverageComparator);        
    }
    
    public void writeStructureVariantV2SortedCoverageReportWithAnnotationToFile(String nameFile , String gffFile, String refFaiFile, int coverageThreshold , int merSize) throws IOException{
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
        String groupCoverageReport = nameFile+".mrkDup.cov.annotation.out";
        
        /**
         * Read annotation file (GFF)
         */
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaiFile, "gene" , merSize);
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        /****************************/

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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
    
    public void writePreciseStructureVariantV2SortedCoverageGroupInfoReportExcel(String nameFile , String refFaiFile, int coverageThreshold) throws IOException{
        /**
        * Suitable for SVGroup object that contain VariaitonV2 object
        * write result to file format for variant report in sorted order (high to low)
        * 
        * "Caution : this function must be called after classifyRoughVariantReport function"
        */
        String tandemReport = nameFile+"_rmDup_tandem.csv";
        String deletionReport = nameFile+"_rmDup_deletion.csv";
        String intraInsertionReport = nameFile+"_rmDup_intraInsertion.csv";
        String interInsertionReport = nameFile+"_rmDup_interInsertion.csv";
        String chimericReport = nameFile+"_rmDup_chimeric.csv";
        String summaryReport = nameFile+"_SummaryStats.txt";
        String groupCoverageReport = nameFile+".mrkDup.cov.PreciseSVType.ginfo.annotation.out";
        
        /**
         * write summary report
         */
        FileWriter writer;
        File f = new File(summaryReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(summaryReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(summaryReport);
        }
        writer.write("Total Breakpoint Found = " + (this.sameChrSVGroup.size()+this.diffChrSVGroup.size())+"\n");
        writer.write("Total BreakPoint with same Chromosome = " + this.sameChrSVGroup.size() + "\n");
        writer.write("Total BreakPoint with different Chromosome = " + this.diffChrSVGroup.size() + "\n");
        writer.write("Possible Tandem = " + this.tandemList.size() + "\n");
        writer.write("Possible Deletion = " + this.deletionList.size()+"\n");
        writer.write("Possible Intra-Insertion = " + this.intraInsertionList.size()+"\n");
        writer.write("Possible Inter-Insertion = " + this.interInsertionList.size()+"\n");
        writer.write("Possible Chimeric = " + this.chimericList.size()+"\n");
        writer.flush();
        writer.close();

        /**
         * Write tandem report
         */
        f = new File(tandemReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(tandemReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(tandemReport);
        }
        
        /**
         * Tandem
         */
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,Correctness(High is best),Annotation Front,Annotation Back"); // header
        writer.write("\n");
        int count = 1;
        for(int i=0;i<this.tandemList.size();i++){
            SVGroup svGroup = this.tandemList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }

            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelTandemSummary());
            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /******************************/
        
        /**
         * Deletion
         */
        
        f = new File(deletionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(deletionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(deletionReport);
        }
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,Correctness(Low is best),Deletion Size,Annotation Front,Annotation Back"); // header
        writer.write("\n");
        count = 1;
        for(int i=0;i<this.deletionList.size();i++){
            SVGroup svGroup = this.deletionList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelDeletionSummary());
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /**********************************/
        
        /**
         * Intra Insertion
         */
        f = new File(intraInsertionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(intraInsertionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(intraInsertionReport);
        }
        writer.write("Group,Chromosome Front(Front Junction),BreakPoint Front(Front Junction),Strand Front(Front Junction)"
                + ",Chromosome Back(Front Junction),BreakPoint Back(Front Junction),Strand Back(Front Junction),Support Reads(Front Junction)"
                + ",Chromosome Front(Back Junction),BreakPoint Front(Back Junction),Strand Front(Back Junction)"
                + ",Chromosome Back(Back Junction),BreakPoint Back(Back Junction),Strand Back(Back Junction),Support Reads(Back Junction)"
                + ",Insertion Size,Insertion Junction,Annotation Front(Front Junction),Annotation Back(Front Junction),Annotation Front(Back Junction),Annotation Back(Back Junction)"); // header
        writer.write("\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.intraTransPairList.entrySet()){
        for(int i=0;i<this.intraInsertionList.size();i++){
            String anno = "";
            SVGroupPair svPair = this.intraInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();
            
            writer.write(count+","+frontJunction.excelShortSummary());
            writer.write(","+backJunction.excelShortSummary());
            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /*****************************************/
        
        /**
         * Inter Translocation
         */
        f = new File(interInsertionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(interInsertionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(interInsertionReport);
        }
        writer.write("Group,Chromosome Front(Front Junction),BreakPoint Front(Front Junction),Strand Front(Front Junction)"
                + ",Chromosome Back(Front Junction),BreakPoint Back(Front Junction),Strand Back(Front Junction),Support Reads(Front Junction)"
                + ",Chromosome Front(Back Junction),BreakPoint Front(Back Junction),Strand Front(Back Junction)"
                + ",Chromosome Back(Back Junction),BreakPoint Back(Back Junction),Strand Back(Back Junction),Support Reads(Back Junction)"
                + ",Insertion Size,Insertion Junction,Annotation Front(Front Junction),Annotation Back(Front Junction),Annotation Front(Back Junction),Annotation Back(Back Junction)"); // header
        writer.write("\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.interTransPairList.entrySet()){
        for(int i=0;i<this.interInsertionList.size();i++){
            String anno = "";
            SVGroupPair svPair = this.interInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();

            writer.write(count+","+frontJunction.excelShortSummary());
            writer.write(","+backJunction.excelShortSummary());
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /**********************************************/
        
        /**
         * Chimeric
         */
        f = new File(chimericReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(chimericReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(chimericReport);
        }
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,Annotation Front,Annotation Back"); // header
        writer.write("\n");
        count = 1;
        Collections.sort(this.chimericList,SVGroup.CoverageComparator);
        for(int i=0;i<this.chimericList.size();i++){
            SVGroup svGroup = this.chimericList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }

            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelShortSummary());
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /***********************************************/
        
        
        /***************************************************/
    }
    
    public void writePreciseStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationExcel(String readFile , String gffFile, String refFaiFile, int coverageThreshold , int merSize) throws IOException{
        /**
        * Suitable for SVGroup object that contain VariaitonV2 object
        * write result to file format for variant report in sorted order (high to low)
        * 
        * "Caution : this function must be called after classifyRoughVariantReport function"
        */
        
        Path path = Paths.get(readFile);
        String fileName = path.getFileName().toString();
        String[] dummy2 = fileName.split("_unmap");         // For more generic this should be split by . (and whole protocal should use . as a filed peaparation like A.unmap.emdup.bam So, we can get read name just split by .)
        String sampleName = dummy2[0];
        
        String tandemReport = path.getParent().toString()+File.separator+"Tandem"+File.separator+sampleName+"_tandem.csv";
        File file = new File(tandemReport);
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParent());
            }
        }
        
        String deletionReport = path.getParent().toString()+File.separator+"Deletion"+File.separator+sampleName+"_deletion.csv";
        file = new File(deletionReport);
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParent());
            }
        }
        
        String intraInsertionReport = path.getParent().toString()+File.separator+"IntraInsertion"+File.separator+sampleName+"_intraInsertion.csv";
        file = new File(intraInsertionReport);
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParent());
            }
        }
        
        String interInsertionReport = path.getParent().toString()+File.separator+"InterInsertion"+File.separator+sampleName+"_interInsertion.csv";
        file = new File(interInsertionReport);
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParent());
            }
        }
        
        String chimericReport = path.getParent().toString()+File.separator+"Chimeric"+File.separator+sampleName+"_chimeric.csv";
        file = new File(chimericReport);
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParent());
            }
        }
//        
//        String tandemReport = nameFile+"_rmDup_tandem.csv";
//        String deletionReport = nameFile+"_rmDup_deletion.csv";
//        String intraInsertionReport = nameFile+"_rmDup_intraInsertion.csv";
//        String interInsertionReport = nameFile+"_rmDup_interInsertion.csv";
//        String chimericReport = nameFile+"_rmDup_chimeric.csv";
        String summaryReport = path.getParent().toString()+File.separator+sampleName+"_SummaryStats.txt";
//        String groupCoverageReport = nameFile+".mrkDup.cov.PreciseSVType.ginfo.annotation.out";
        
        
        /**
         * write summary report
         */
        FileWriter writer;
        File f = new File(summaryReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(summaryReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(summaryReport);
        }
        writer.write("Total Breakpoint Found = " + (this.sameChrSVGroup.size()+this.diffChrSVGroup.size())+"\n");
        writer.write("Total BreakPoint with same Chromosome = " + this.sameChrSVGroup.size() + "\n");
        writer.write("Total BreakPoint with different Chromosome = " + this.diffChrSVGroup.size() + "\n");
        writer.write("Possible Tandem = " + this.tandemList.size() + "\n");
        writer.write("Possible Deletion = " + this.deletionList.size()+"\n");
        writer.write("Possible Intra-Insertion = " + this.intraInsertionList.size()+"\n");
        writer.write("Possible Inter-Insertion = " + this.interInsertionList.size()+"\n");
        writer.write("Possible Chimeric = " + this.chimericList.size()+"\n");
        writer.flush();
        writer.close();
        
        
        /**
         * Read annotation file (GFF)
         */
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV4(gffFile,refFaiFile, "gene" , merSize);
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        
        System.out.println("refAnno is empty? = "+refAnnoIndex.size());
        /****************************/

        /**
         * Write tandem report
         */
        f = new File(tandemReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(tandemReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(tandemReport);
        }
        
        /**
         * Tandem
         */
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,Correctness(High is best),numStrandPattern,+|+,-|-,+|-,-|+,Annotation Front,Annotation Back"); // header
        writer.write("\n");
        int count = 1;
        for(int i=0;i<this.tandemList.size();i++){
            SVGroup svGroup = this.tandemList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
//            /**
//             * annotation part
//             * Find annotation for this svGroup
//             */
//            
//            if(!svGroup.isIdentityFlag()){
//                svGroup.defineIdentity();
//            }
//            
//            long avgIniPosF = svGroup.getRPF();
//            long avgLastPosB = svGroup.getRPB();
//            
//            long strandChrF = ((long)svGroup.getStrandF()<<11)+(long)svGroup.getNumChrF();
//            long strandChrB = ((long)svGroup.getStrandB()<<11)+(long)svGroup.getNumChrB();
//            
//            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            String strAnnoF = "null";
//            String strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                svGroup.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                svGroup.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
            
            /**
            * Annotation part
            */
            ArrayList<String> annotation = getAnnotationOfSVGroupV3(svGroup, refAnno, refAnnoIndex, "gene",merSize);
            String strAnnoF = annotation.get(0);
            String strAnnoB = annotation.get(1);
            /******************/
                
            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelTandemSummary());
            writer.write(","+strAnnoF+","+strAnnoB);
            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /******************************/
        
        /**
         * Deletion
         */
        
        f = new File(deletionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(deletionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(deletionReport);
        }
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,Correctness(Low is best),Deletion Size,numStrandPattern,+|+,-|-,+|-,-|+,Annotation Front,Annotation Back"); // header
        writer.write("\n");
        count = 1;
        for(int i=0;i<this.deletionList.size();i++){
            SVGroup svGroup = this.deletionList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
//            /**
//             * annotation part
//             * Find annotation for this svGroup
//             */
//            
//            if(!svGroup.isIdentityFlag()){
//                svGroup.defineIdentity();
//            }
//            
//            long avgIniPosF = svGroup.getRPF();
//            long avgLastPosB = svGroup.getRPB();
//            
//            long strandChrF = ((long)svGroup.getStrandF()<<11)+(long)svGroup.getNumChrF();
//            long strandChrB = ((long)svGroup.getStrandB()<<11)+(long)svGroup.getNumChrB();
//            
//            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            String strAnnoF = "null";
//            String strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                svGroup.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                svGroup.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }

            /**
            * Annotation part
            */
            ArrayList<String> annotation = getAnnotationOfSVGroupV3(svGroup, refAnno, refAnnoIndex, "gene",merSize);
            String strAnnoF = annotation.get(0);
            String strAnnoB = annotation.get(1);
            /******************/
            
            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelDeletionSummary());
            writer.write("," + strAnnoF + "," + strAnnoB);
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /**********************************/
        
        /**
         * Intra Insertion
         */
        f = new File(intraInsertionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(intraInsertionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(intraInsertionReport);
        }
        writer.write("Group,Chromosome Front(Front Junction),BreakPoint Front(Front Junction),Strand Front(Front Junction)"
                + ",Chromosome Back(Front Junction),BreakPoint Back(Front Junction),Strand Back(Front Junction),Support Reads(Front Junction),numStrandPattern,+|+,-|-,+|-,-|+"
                + ",Chromosome Front(Back Junction),BreakPoint Front(Back Junction),Strand Front(Back Junction)"
                + ",Chromosome Back(Back Junction),BreakPoint Back(Back Junction),Strand Back(Back Junction),Support Reads(Back Junction),numStrandPattern,+|+,-|-,+|-,-|+"
                + ",Insert Portion Size,Insertion Junction,Annotation Front(Front Junction),Annotation Back(Front Junction),Annotation Front(Back Junction),Annotation Back(Back Junction)"); // header
        writer.write("\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.intraTransPairList.entrySet()){
        for(int i=0;i<this.intraInsertionList.size();i++){
            String anno = "";
            SVGroupPair svPair = this.intraInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();
        
//            ArrayList<SVGroup> intraTransPair = map.getKey();
//            SVGroup frontJunction = intraTransPair.get(0);
//            SVGroup backJunction = intraTransPair.get(1);
            
//             /**
//             * annotation part Front Junction
//             * Find annotation for this svGroup
//             */
//            
//            if(!frontJunction.isIdentityFlag()){
//                frontJunction.defineIdentity();
//            }
//            
//            long avgIniPosF = frontJunction.getRPF();
//            long avgLastPosB = frontJunction.getRPB();
//
//            long strandChrF = ((long)frontJunction.getStrandF()<<11)+(long)frontJunction.getNumChrF();
//            long strandChrB = ((long)frontJunction.getStrandB()<<11)+(long)frontJunction.getNumChrB();
//            
//            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopF = ((strandChrF<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            long chrPosStartB = ((strandChrB<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            long chrPosStartF = (((long)frontJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopF = (((long)frontJunction.getNumChrF()<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            long chrPosStartB = (((long)frontJunction.getNumChrB()<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopB = (((long)frontJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            String strAnnoF = "null";
//            String strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                frontJunction.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                frontJunction.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
//            /**********************/
            /**
            * Annotation part
            */
            ArrayList<String> annotationF = getAnnotationOfSVGroupV3(frontJunction, refAnno, refAnnoIndex, "gene", merSize);
            String strAnnoF = annotationF.get(0);
            String strAnnoB = annotationF.get(1);
            /******************/
            

            writer.write(count+","+frontJunction.excelShortSummary());
            anno = anno+strAnnoF+","+strAnnoB;
//            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
//            writer.write("\n");
            
//            /**
//             * annotation part back Junction
//             * Find annotation for this svGroup
//             */
            
//            if(!backJunction.isIdentityFlag()){
//                backJunction.defineIdentity();
//            }
//            
//            avgIniPosF = backJunction.getRPF();
//            avgLastPosB = backJunction.getRPB();
//            
//            strandChrF = ((long)backJunction.getStrandF()<<11)+(long)backJunction.getNumChrF();
//            strandChrB = ((long)backJunction.getStrandB()<<11)+(long)backJunction.getNumChrB();
//            
//            chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            chrPosStopF = ((strandChrF<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            chrPosStartB = ((strandChrB<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            chrPosStartF = (((long)backJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            chrPosStopF = (((long)backJunction.getNumChrF()<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            chrPosStartB = (((long)backJunction.getNumChrB()<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            chrPosStopB = (((long)backJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            strAnnoF = "null";
//            strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                backJunction.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                backJunction.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
//            /**********************/
            /**
            * Annotation part
            */
            ArrayList<String> annotationB = getAnnotationOfSVGroupV3(backJunction, refAnno, refAnnoIndex, "gene", merSize);
            strAnnoF = annotationB.get(0);
            strAnnoB = annotationB.get(1);
            /******************/
            
            writer.write(","+backJunction.excelShortSummary());
            anno = anno + "," + strAnnoF+","+strAnnoB;
            writer.write("," + svPair.getInsertionSize()+","+svPair.getInsertionJunction()+","+ anno);
            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /*****************************************/
        
        /**
         * Inter Translocation
         */
        f = new File(interInsertionReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(interInsertionReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(interInsertionReport);
        }
        writer.write("Group,Chromosome Front(Front Junction),BreakPoint Front(Front Junction),Strand Front(Front Junction)"
                + ",Chromosome Back(Front Junction),BreakPoint Back(Front Junction),Strand Back(Front Junction),Support Reads(Front Junction),numStrandPattern,+|+,-|-,+|-,-|+"
                + ",Chromosome Front(Back Junction),BreakPoint Front(Back Junction),Strand Front(Back Junction)"
                + ",Chromosome Back(Back Junction),BreakPoint Back(Back Junction),Strand Back(Back Junction),Support Reads(Back Junction),numStrandPattern,+|+,-|-,+|-,-|+"
                + ",Insert Portion Size,Insertion Junction,Annotation Front(Front Junction),Annotation Back(Front Junction),Annotation Front(Back Junction),Annotation Back(Back Junction)"); // header
        writer.write("\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.interTransPairList.entrySet()){
        for(int i=0;i<this.interInsertionList.size();i++){
            String anno = "";
            SVGroupPair svPair = this.interInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();
            
//            ArrayList<SVGroup> intraTransPair = map.getKey();
//            SVGroup frontJunction = intraTransPair.get(0);
//            SVGroup backJunction = intraTransPair.get(1);
            
//             /**
//             * annotation part Front Junction
//             * Find annotation for this svGroup
//             */
//            
//            if(!frontJunction.isIdentityFlag()){
//                frontJunction.defineIdentity();
//            }
//            
//            long avgIniPosF = frontJunction.getRPF();
//            long avgLastPosB = frontJunction.getRPB();
//
//            long strandChrF = ((long)frontJunction.getStrandF()<<11)+(long)frontJunction.getNumChrF();
//            long strandChrB = ((long)frontJunction.getStrandB()<<11)+(long)frontJunction.getNumChrB();
//            
//            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopF = ((strandChrF<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            long chrPosStartB = ((strandChrB<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            long chrPosStartF = (((long)frontJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopF = (((long)frontJunction.getNumChrF()<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            long chrPosStartB = (((long)frontJunction.getNumChrB()<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopB = (((long)frontJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            String strAnnoF = "null";
//            String strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                frontJunction.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                frontJunction.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
//            /**********************/
            /**
            * Annotation part
            */
            ArrayList<String> annotationF = getAnnotationOfSVGroupV3(frontJunction, refAnno, refAnnoIndex, "gene", merSize);
            String strAnnoF = annotationF.get(0);
            String strAnnoB = annotationF.get(1);
            /******************/
            
            writer.write(count+","+frontJunction.excelShortSummary());
            anno = anno + strAnnoF+","+strAnnoB;
//            writer.write("\t" + anno);
//            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);

            
//            /**
//             * annotation part back Junction
//             * Find annotation for this svGroup
//             */
//            
//            if(!backJunction.isIdentityFlag()){
//                backJunction.defineIdentity();
//            }
//            
//            avgIniPosF = backJunction.getRPF();
//            avgLastPosB = backJunction.getRPB();
//            
//            strandChrF = ((long)backJunction.getStrandF()<<11)+(long)backJunction.getNumChrF();
//            strandChrB = ((long)backJunction.getStrandB()<<11)+(long)backJunction.getNumChrB();
//            
//            chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            chrPosStopF = ((strandChrF<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            chrPosStartB = ((strandChrB<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            chrPosStartF = (((long)backJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            chrPosStopF = (((long)backJunction.getNumChrF()<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            chrPosStartB = (((long)backJunction.getNumChrB()<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            chrPosStopB = (((long)backJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            strAnnoF = "null";
//            strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                backJunction.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                backJunction.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
//            /**********************/
            /**
            * Annotation part
            */
            ArrayList<String> annotationB = getAnnotationOfSVGroupV3(backJunction, refAnno, refAnnoIndex, "gene", merSize);
            strAnnoF = annotationB.get(0);
            strAnnoB = annotationB.get(1);
            /******************/

            writer.write(","+backJunction.excelShortSummary());
            anno = anno + "," + strAnnoF+","+strAnnoB;
            writer.write(","+svPair.getInsertionSize()+","+svPair.getInsertionJunction()+ "," + anno);
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /**********************************************/
        
        /**
         * Chimeric
         */
        f = new File(chimericReport); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(chimericReport,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(chimericReport);
        }
        writer.write("Group,Chromosome Front,BreakPoint Front,Strand Front,Chromosome Back,BreakPoint Back,Strand Back,Support Reads,numStrandPattern,+|+,-|-,+|-,-|+,Annotation Front,Annotation Back"); // header
        writer.write("\n");
        count = 1;
        Collections.sort(this.chimericList,SVGroup.CoverageComparator);
        for(int i=0;i<this.chimericList.size();i++){
            SVGroup svGroup = this.chimericList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
//            /**
//             * annotation part
//             * Find annotation for this svGroup
//             */
//            
//            if(!svGroup.isIdentityFlag()){
//                svGroup.defineIdentity();
//            }
//            
//            long avgIniPosF = svGroup.getRPF();
//            long avgLastPosB = svGroup.getRPB();
//
//            long strandChrF = ((long)svGroup.getStrandF()<<11)+(long)svGroup.getNumChrF();
//            long strandChrB = ((long)svGroup.getStrandB()<<11)+(long)svGroup.getNumChrB();
//            
//            long chrPosStartF = ((strandChrF<<28)+avgIniPosF)<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopF = ((strandChrF<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            long chrPosStartB = ((strandChrB<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [strand 1bit][chr11bit][position28bit][empty23bit] 
//            long chrPosStopB = ((strandChrB<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//            
////            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
////            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
////            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
//
//            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
//            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);
//
//            String strAnnoF = "null";
//            String strAnnoB = "null";
//            if(annoGroupIndexF >= 0){
//                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
//                svGroup.setAnnoF(annoF);
//                strAnnoF = annoF.toString();
//            }
//            if(annoGroupIndexB >= 0){
//                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
//                svGroup.setAnnoB(annoB);
//                strAnnoB = annoB.toString();
//            }
//
//            /**********************/
            
            /**
            * Annotation part
            */
            ArrayList<String> annotation = getAnnotationOfSVGroupV3(svGroup, refAnno, refAnnoIndex, "gene",merSize);
            String strAnnoF = annotation.get(0);
            String strAnnoB = annotation.get(1);
            /******************/

            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(count+","+svGroup.excelShortSummary());
            writer.write("," + strAnnoF + "," + strAnnoB);
            writer.write("\n");
            count++;
        }
        writer.flush();
        writer.close();
        /***********************************************/
        
        
        /***************************************************/
    }
    
    public void writePreciseStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(String nameFile , String gffFile, String refFaiFile, int coverageThreshold , int merSize) throws IOException{
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
        String groupCoverageReport = nameFile+".mrkDup.cov.PreciseSVType.ginfo.annotation.out";
        
        /**
         * Read annotation file (GFF)
         */
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaiFile, "gene" , merSize);
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        /****************************/

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
        
        /**
         * Tandem
         */
        writer.write("Tandem\n");
        int count = 1;
        for(int i=0;i<this.tandemList.size();i++){
            SVGroup svGroup = this.tandemList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.shortTandemSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            count++;
        }
        
        /******************************/
        
        /**
         * Deletion
         */
        writer.write("Deletion\n");
        count = 1;
        for(int i=0;i<this.deletionList.size();i++){
            SVGroup svGroup = this.deletionList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.shortDeletionSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            count++;
        }
        /**********************************/
        
        /**
         * Intra Insertion
         */
        writer.write("Intra-Insertion\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.intraTransPairList.entrySet()){
        for(int i=0;i<this.intraInsertionList.size();i++){
            SVGroupPair svPair = this.intraInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();
        
//            ArrayList<SVGroup> intraTransPair = map.getKey();
//            SVGroup frontJunction = intraTransPair.get(0);
//            SVGroup backJunction = intraTransPair.get(1);
            
             /**
             * annotation part Front Junction
             * Find annotation for this svGroup
             */
            
            if(!frontJunction.isIdentityFlag()){
                frontJunction.defineIdentity();
            }
            
            long avgIniPosF = frontJunction.getRPF();
            long avgLastPosB = frontJunction.getRPB();

            long chrPosStartF = (((long)frontJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)frontJunction.getNumChrF()<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)frontJunction.getNumChrB()<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)frontJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                frontJunction.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                frontJunction.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }
            /**********************/

            writer.write(">"+count+"\t"+frontJunction.shortSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            
            /**
             * annotation part back Junction
             * Find annotation for this svGroup
             */
            
            if(!backJunction.isIdentityFlag()){
                backJunction.defineIdentity();
            }
            
            avgIniPosF = backJunction.getRPF();
            avgLastPosB = backJunction.getRPB();

            chrPosStartF = (((long)backJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            chrPosStopF = (((long)backJunction.getNumChrF()<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            chrPosStartB = (((long)backJunction.getNumChrB()<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            chrPosStopB = (((long)backJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            strAnnoF = "null";
            strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                backJunction.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                backJunction.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }
            /**********************/
            writer.write("\t"+backJunction.shortSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            count++;
        }             
        /*****************************************/
        
        /**
         * Inter Translocation
         */
        writer.write("Inter-Insertion\n");
        count = 1;
//        for(Map.Entry<ArrayList<SVGroup>,Boolean> map : this.interTransPairList.entrySet()){
        for(int i=0;i<this.interInsertionList.size();i++){
            SVGroupPair svPair = this.interInsertionList.get(i);
            SVGroup frontJunction = svPair.getFrontSVGroup();
            SVGroup backJunction = svPair.getBackSVGroup();
            
//            ArrayList<SVGroup> intraTransPair = map.getKey();
//            SVGroup frontJunction = intraTransPair.get(0);
//            SVGroup backJunction = intraTransPair.get(1);
            
             /**
             * annotation part Front Junction
             * Find annotation for this svGroup
             */
            
            if(!frontJunction.isIdentityFlag()){
                frontJunction.defineIdentity();
            }
            
            long avgIniPosF = frontJunction.getRPF();
            long avgLastPosB = frontJunction.getRPB();

            long chrPosStartF = (((long)frontJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)frontJunction.getNumChrF()<<28)+frontJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)frontJunction.getNumChrB()<<28)+frontJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)frontJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                frontJunction.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                frontJunction.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }
            /**********************/

            writer.write(">"+count+"\t"+frontJunction.shortSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            
            /**
             * annotation part back Junction
             * Find annotation for this svGroup
             */
            
            if(!backJunction.isIdentityFlag()){
                backJunction.defineIdentity();
            }
            
            avgIniPosF = backJunction.getRPF();
            avgLastPosB = backJunction.getRPB();

            chrPosStartF = (((long)backJunction.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            chrPosStopF = (((long)backJunction.getNumChrF()<<28)+backJunction.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            chrPosStartB = (((long)backJunction.getNumChrB()<<28)+backJunction.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            chrPosStopB = (((long)backJunction.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            strAnnoF = "null";
            strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                backJunction.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                backJunction.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }
            /**********************/
            writer.write("\t"+backJunction.shortSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            count++;
        }             
        /**********************************************/
        
        /**
         * Inter Translocation
         */
        writer.write("Chimeric\n");
        count = 1;
        Collections.sort(this.chimericList,SVGroup.CoverageComparator);
        for(int i=0;i<this.chimericList.size();i++){
            SVGroup svGroup = this.chimericList.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                continue;
            }
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.shortSummary());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");
            count++;
        }
        /***********************************************/
        
        writer.flush();
        writer.close();
        /***************************************************/
    }
    
    
    public void writeStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(String nameFile , String gffFile, String refFaiFile, int coverageThreshold , int merSize) throws IOException{
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
        String groupCoverageReport = nameFile+".mrkDup.cov.ginfo.annotation.out";
        
        /**
         * Read annotation file (GFF)
         */
        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaiFile, "gene" , merSize);
        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        /****************************/

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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
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
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            long avgIniPosF = svGroup.getRPF();
            long avgLastPosB = svGroup.getRPB();

            long chrPosStartF = (((long)svGroup.getNumChrF()<<28)+avgIniPosF)<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopF = (((long)svGroup.getNumChrF()<<28)+svGroup.getRPF())<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)
            long chrPosStartB = (((long)svGroup.getNumChrB()<<28)+svGroup.getRPB())<<23;     // create chrPosStart code [chr5bit][position28bit][empty23bit] 
            long chrPosStopB = (((long)svGroup.getNumChrB()<<28)+avgLastPosB)<<23;     // (the reason that we have to have empty bit 23 bit at the back is It help us to do binary search more easily with reference annotation. the reference annotation have been operate by AND with mask that will make the 23bit on the back chenge to 0 value)

            int annoGroupIndexF = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartF, chrPosStopF);
            int annoGroupIndexB = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStartB, chrPosStopB);

            String strAnnoF = "null";
            String strAnnoB = "null";
            if(annoGroupIndexF >= 0){
                Annotation annoF = refAnnoIndex.get(annoGroupIndexF);
                svGroup.setAnnoF(annoF);
                strAnnoF = annoF.toString();
            }
            if(annoGroupIndexB >= 0){
                Annotation annoB = refAnnoIndex.get(annoGroupIndexB);
                svGroup.setAnnoB(annoB);
                strAnnoB = annoB.toString();
            }

            /**********************/
            
            ArrayList<VariationV2> varList = svGroup.getVarList();
            writer.write(">"+count+"\t"+svGroup.toString());
            writer.write("\t" + "Annotation Front : " + strAnnoF + "\t" + "Annotation back : " + strAnnoB);
            writer.write("\n");            
            count++;
        }
        
        writer.flush();
        writer.close();
        /***************************************************/
    }
    
    public void writeStructureVariantV2EuclideanDistanceTable(String nameFile, int coverageThreshold) throws IOException{
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
        String groupCoverageReport = nameFile+".mrkDup.cov.ginfo.annotation.out";
        String euclideanFile = nameFile+"euclideanDistance.out";
        
        /**
         * Read annotation file (GFF)
         */
//        ReferenceAnnotation refAnno = SequenceUtil.readAnnotationFileV3(gffFile,refFaiFile, "gene" , merSize);
//        Map<Integer,Annotation> refAnnoIndex = refAnno.getAnnotationIndex();
        /****************************/

        /**
         * Write tandem cov report
         */
        FileWriter writer;
        File f = new File(euclideanFile); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(euclideanFile,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(euclideanFile);
        }
        
        writer.write("Same Chromosome Euclidean distance table\n");
        int count = 1;
        for(int i=0;i<this.sameChrSVGroup.size();i++){
            SVGroup svGroup = this.sameChrSVGroup.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                break;
            }
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            writer.write(">"+count+"\t"+svGroup.shortSummary()+"\t|");
            ArrayList<Double> euList = svGroup.getSameChrEnclideanDistance();
            /**
             * print All euclidean value
             */            
//            for(int k=0;k<euList.size();k++){
//                DecimalFormat df = new DecimalFormat("#.00"); // limit double to 2 decimal
//                writer.write("\t"+df.format(euList.get(k)));
//                int minIndex = euList.indexOf(Collections.min(euList));
//            }
            /*************/
            
            /**
             * pick top three lowest euclidean distance value
             */
            ArrayList<Double> dummyEuList = new ArrayList<Double>(euList);
            for(int k=0;k<3;k++){
                int minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                DecimalFormat df = new DecimalFormat("#.00"); // limit double to 2 decimal
                double minEuDis = dummyEuList.get(minIndex);
                int realMinIndex = euList.indexOf(minEuDis);        // true min index for write report
                if(realMinIndex+1 != count){ // not report eu distance that calculate between it selfe 
                    writer.write("\t["+(realMinIndex+1)+"]"+df.format(dummyEuList.get(minIndex))); // realminindex plus one because index when write is start at 1 but real index is start at 0
                }
                dummyEuList.remove(minIndex);
            }
            /************/
            
            writer.write("\n");
            count++;
        }
        
        writer.write("Difference Chromosome Euclidean distance table\n");
        count = 1;
        for(int i=0;i<this.diffChrSVGroup.size();i++){
            SVGroup svGroup = this.diffChrSVGroup.get(i);
            if(svGroup.getNumCoverage()<coverageThreshold){
                // fall threshold ignore this group
                break;
            }
            
            /**
             * annotation part
             * Find annotation for this svGroup
             */
            
            if(!svGroup.isIdentityFlag()){
                svGroup.defineIdentity();
            }
            
            writer.write(">"+count+"\t"+svGroup.shortSummary()+"\t|");
            ArrayList<Double> euList = svGroup.getDiffChrEuclideanDistance();
            
            /**
             * print All euclidean value
             */ 
//            for(int k=0;k<euList.size();k++){
//                DecimalFormat df = new DecimalFormat("#.00"); // limit double to 2 decimal
//                writer.write("\t"+df.format(euList.get(k)));
//            }
            /**********/
            
            /**
             * pick top three lowest euclidean distance value
             */
            ArrayList<Double> dummyEuList = new ArrayList<Double>(euList);
            for(int k=0;k<3;k++){
                int minIndex = dummyEuList.indexOf(Collections.min(dummyEuList));
                DecimalFormat df = new DecimalFormat("#.00"); // limit double to 2 decimal
                double minEuDis = dummyEuList.get(minIndex);
                int realMinIndex = euList.indexOf(minEuDis);    // true min index for write report
                if(realMinIndex+1 != count){
                    writer.write("\t["+(realMinIndex+1)+"]"+df.format(dummyEuList.get(minIndex)));  // realminindex plus one because index when write is start at 1 but real index is start at 0
                }
                dummyEuList.remove(minIndex);
            }
            /************/

            writer.write("\n");
            count++;
        }
        
        writer.flush();
        writer.close();
        /***************************************************/
    }
    
    public void createRefIndex(String refFaIdxFile) throws FileNotFoundException, IOException{
        RandomAccessFile rbRefIdx = new RandomAccessFile(refFaIdxFile,"r");
        String line;
        String name = "";
        long numChr = 0;
        
        while ((line = rbRefIdx.readLine()) != null) {
            String[] linePortion = line.split("\t");
            name = linePortion[0];
            this.refIndex.put(name, ++numChr);
        }
        rbRefIdx.close();
    }
    
    public String identifySameChrPreciseSVType(SVGroup main,SVGroup sub){
        if(main.getStrandF()==0 && main.getStrandB()==0){
            // main has ++ strand
            if(sub.getStrandF()==0 && sub.getStrandB()==0){
//                if(main.getRPF() == 140773609 && main.getRPB() == 140777194){
//                    System.out.println();
//                }else if(main.getRPF() == 140772667 && main.getRPB() == 140773504){
//                    System.out.println();
//                }
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();                            
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.intraInsertionList.add(svPair);
                        }
                        // classify main as tandem or delete in the same time when class as intrainsertion 
                        if(main.getRPF()>main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
//                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
//                            return "deletion";
                        }
                        
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.intraInsertionList.add(svPair);
                        }
                        
                        // classify main as tandem or delete in the same time when class as intrainsertion 
                        if(main.getRPF()>main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
//                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
//                            return "deletion";
                        }
                        
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else{
                        if(main.getRPF()>main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
                            return "deletion";
                        }
                    }
                }else{
                    if(main.getRPF()>main.getRPB()){
                        this.tandemList.add(main);
                        main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
                        return "tandem";
                    }else{
                        this.deletionList.add(main);
                        main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
                        return "deletion";
                    }
                }
            }else{
                if(main.getRPF()>main.getRPB()){
                    this.tandemList.add(main);
                    main.setTandemFlag(true);
//                    main.setSvType("tandem");
//                    main.setSvTypeCode((byte)0);
                    return "tandem";
                }else{
                    this.deletionList.add(main);
                    main.setDeleteFlag(true);
//                    main.setSvType("deletion");
//                    main.setSvTypeCode((byte)1);
                    return "deletion";
                }
            }
        }else if(main.getStrandF()==1 && main.getStrandB()==1){
            // main has -- strand
            if(sub.getStrandF()==1 && sub.getStrandB()==1){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.intraInsertionList.add(svPair);
                        }
                        
                        // classify main as tandem or delete in the same time when class as intrainsertion 
                        if(main.getRPF()<main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
//                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
//                            return "deletion";
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.intraInsertionList.add(svPair);
                        }
                        
                        // classify main as tandem or delete in the same time when class as intrainsertion 
                        if(main.getRPF()<main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
//                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
//                            return "deletion";
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else{
                        if(main.getRPF()<main.getRPB()){
                            this.tandemList.add(main);
                            main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
                            return "tandem";
                        }else{
                            this.deletionList.add(main);
                            main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
                            return "deletion";
                        }
                    }
                }else{
                    if(main.getRPF()<main.getRPB()){
                        this.tandemList.add(main);
                        main.setTandemFlag(true);
//                            main.setSvType("tandem");
//                            main.setSvTypeCode((byte)0);
                        return "tandem";
                    }else{
                        this.deletionList.add(main);
                        main.setDeleteFlag(true);
//                            main.setSvType("deletion");
//                            main.setSvTypeCode((byte)1);
                        return "deletion";
                    }
                }
            }else{
                if(main.getRPF()<main.getRPB()){
                    this.tandemList.add(main);
                    main.setTandemFlag(true);
//                    main.setSvType("tandem");
//                    main.setSvTypeCode((byte)0);
                    return "tandem";
                }else{
                    this.deletionList.add(main);
                    main.setDeleteFlag(true);
//                    main.setSvType("deletion");
//                    main.setSvTypeCode((byte)1);
                    return "deletion";
                }
            }
        }else if(main.getStrandF()==0 && main.getStrandB()==1){
            if(sub.getStrandF()==1 && sub.getStrandB()==0){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.intraInsertionList.add(svPair);
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.intraInsertionList.add(svPair);
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else{
                        // try with mext min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }else if(main.getStrandF()==1 && main.getStrandB()==0){
            if(sub.getStrandF()==0 && sub.getStrandB()==1){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.intraInsertionList.add(svPair);
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.intraTransPairList.put(dummyList, true);
//                        main.setIntraTransFlag(true);
//                        sub.setIntraTransFlag(true);
                        if(!(main.isIntraInsertionFlag() && sub.isIntraInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setIntraInsertionFlag(true);
                            sub.setIntraInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.intraInsertionList.add(svPair);
                        }
//                        main.setSvType("intraTrans");
//                        main.setSvTypeCode((byte)2);
//                        sub.setSvType("intraTrans");
//                        sub.setSvTypeCode((byte)2);
                        return "intraIns";
                    }else{
                        // try with mext min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }        
        return null;
    }
    
    public String identifyDiffChrPreciseSVType(SVGroup main,SVGroup sub){       
        if(main.getStrandF()==0 && main.getStrandB()==0){
            // main has ++ strand
            if(sub.getStrandF()==0 && sub.getStrandB()==0){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.interInsertionList.add(svPair);
                        }
                        
//                        this.interTransPairList.put(dummyList, true);
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.interInsertionList.add(svPair);
                        }
//                        this.interTransPairList.put(dummyList, true);                        
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else{
                        // try with next min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
//                    main.setSvType("chimeric");
//                    main.setSvTypeCode((byte)4);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }else if(main.getStrandF()==1 && main.getStrandB()==1){
            // main has -- strand
            if(sub.getStrandF()==1 && sub.getStrandB()==1){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.interInsertionList.add(svPair);
                        }
                        
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)    
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.interInsertionList.add(svPair);
                        }
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else{
                        // try with next min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
//                    main.setSvType("chimeric");
//                    main.setSvTypeCode((byte)4);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }else if(main.getStrandF()==0 && main.getStrandB()==1){
            if(sub.getStrandF()==1 && sub.getStrandB()==0){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.interInsertionList.add(svPair);
                        }
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)    
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.interInsertionList.add(svPair);
                        }
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else{
                        // try with next min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
//                    main.setSvType("chimeric");
//                    main.setSvTypeCode((byte)4);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }else if(main.getStrandF()==1 && main.getStrandB()==0){
            if(sub.getStrandF()==0 && sub.getStrandB()==1){
                if(main.getChrF().equals(sub.getChrB()) && main.getChrB().equals(sub.getChrF())){
                    if(main.getRPF()<sub.getRPB() && main.getRPB()>sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(sub);
//                        dummyList.add(main);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(sub);
                            svPair.setBackSVGroup(main);
                            this.interInsertionList.add(svPair);
                        }
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else if(main.getRPF()>sub.getRPB() && main.getRPB()<sub.getRPF()){
//                        ArrayList<SVGroup> dummyList = new ArrayList();
//                        dummyList.add(main);
//                        dummyList.add(sub);
//                        this.interTransPairList.put(dummyList, true);
//                        main.setInterTransFlag(true);
//                        sub.setInterTransFlag(true);
                        if(!(main.isInterInsertionFlag() && sub.isInterInsertionFlag())){
                            // add to list (only svgroup that did not have add before)
                            SVGroupPair svPair = new SVGroupPair();
                            main.setInterInsertionFlag(true);
                            sub.setInterInsertionFlag(true);
                            svPair.setFrontSVGroup(main);
                            svPair.setBackSVGroup(sub);
                            this.interInsertionList.add(svPair);
                        }
//                        main.setSvType("interTrans");
//                        main.setSvTypeCode((byte)3);
//                        sub.setSvType("interTrans");
//                        sub.setSvTypeCode((byte)3);
                        return "interIns";
                    }else{
                        // try with next min
                        return null;
                    }
                }else{
                    this.chimericList.add(main);
                    main.setChimericFlag(true);
//                    main.setSvType("chimeric");
//                    main.setSvTypeCode((byte)4);
                    return "chimeric";
                }
            }else{
                this.chimericList.add(main);
                main.setChimericFlag(true);
//                main.setSvType("chimeric");
//                main.setSvTypeCode((byte)4);
                return "chimeric";
            }
        }
        return null;
    }
    
    public void identifyCorrectness(String bamFile,String samtoolsDirectory){
        /**
         * Analyze correctness of SVType (tandem and delete) for each svGroup 
         * 
         * Loop each svType list and check with our ideal hypothesis we call "relation of coverage with correctness"
         * 
         * Utilize samtools view command to get coverage of each specific region (our breakpoint)
         */
        
        String baseCommand = samtoolsDirectory + " view -F 0X100 " + bamFile;                // base command is samtools view with option 0F 0x100 which mean not include secondary map read
        
        // Loop deletion list
        // the correctness of deletion is the coverage of base that has position in the middle of deletion (the coverage should be zero if it is real deletion)
        for(int i=0;i<this.deletionList.size();i++){
            String addCommand = "";
            SVGroup svGroup = deletionList.get(i);
            
            // query Coverage of middle base it the deletiob region
            int deletionSize = svGroup.getDeletionSize();
            long breakpointF = svGroup.getRPF();
            long middleDelBasePos = breakpointF + (deletionSize/2);
            String chrNameF = svGroup.getChrF();
            addCommand = "chr"+chrNameF+":"+middleDelBasePos+"-"+middleDelBasePos;
            String command = baseCommand + " " + addCommand + " | wc";
            String result = executeCommandLine(command);
            String[] splitRes = result.split("\\s+");
            int coverage = Integer.parseInt(splitRes[1]);
            svGroup.setDeletionCorrectness(coverage);
            
//            // query Coverage of frontbreakpoint
//            long breakpointF = svGroup.getRPF();
//            String chrNameF = svGroup.getChrF();
//            addCommand = "chr"+chrNameF+":"+breakpointF+"-"+breakpointF;
//            String command = baseCommand + " " + addCommand + " | wc";
            
//            String resultF = executeCommandLine(command);
//            String[] splitResF = resultF.split("\\s+");
//            int coverageF = Integer.parseInt(splitResF[1]);
//            
//            // query Coverage f backreakpoint
//            long breakpointB = svGroup.getRPB();
//            String chrNameB = svGroup.getChrF();
//            addCommand = "chr"+chrNameB+":"+breakpointB+"-"+breakpointB;
//            command = baseCommand + " " + addCommand + " | wc";
//            
//            String resultB = executeCommandLine(command);
//            String[] splitResB = resultB.split("\\s+");
//            int coverageB = Integer.parseInt(splitResB[1]);
//            
//            // Add sum of two coverage to delectioncorrectness field if SVGroup
//            svGroup.setDeletionCorrectness(coverageF + coverageB);   
        }
        
        
        // Loop tandem List
        for(int i=0;i<this.tandemList.size();i++){
            String addCommand = "";
            SVGroup svGroup = tandemList.get(i);
            
            // query Coverage of frontbreakpoint
            long breakpointF = svGroup.getRPF();
            String chrNameF = svGroup.getChrF();
            addCommand = "chr"+chrNameF+":"+breakpointF+"-"+breakpointF;
            String command = baseCommand + " " + addCommand + " | wc";
            
            String resultF = executeCommandLine(command);
            String[] splitResF = resultF.split("\\s+");
            int coverageF = Integer.parseInt(splitResF[1]);
            
            // query Coverage f backreakpoint
            long breakpointB = svGroup.getRPB();
            String chrNameB = svGroup.getChrF();
            addCommand = "chr"+chrNameB+":"+breakpointB+"-"+breakpointB;
            command = baseCommand + " " + addCommand + " | wc";
            
            String resultB = executeCommandLine(command);
            String[] splitResB = resultB.split("\\s+");
            int coverageB = Integer.parseInt(splitResB[1]);
            
            // Add sum of two coverage to tandemcorrectness field if SVGroup
            svGroup.setTandemCorrectness(coverageF + coverageB);   
        }

        Collections.sort(this.deletionList,SVGroup.CoverageCorrectnessDeletionComparator);
        Collections.sort(this.tandemList,SVGroup.CoverageCorrectnessTandemComparator);
        /********************/
    }
    
    private String executeCommandLine(String command) {

            StringBuffer output = new StringBuffer();
            String[] cmd = {"/bin/sh","-c",command};
            Process p;
            try {
                    p = Runtime.getRuntime().exec(cmd);
                    p.waitFor();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));

                    String line = "";			
                    while ((line = reader.readLine())!= null) {
                        output.append(line + "\n");
                    }

            } catch (Exception e) {
                    e.printStackTrace();
            }

            return output.toString();
    }
    
    public String getChrNameFromRefIndex(Map<String, Long> refIndexMap, Long chrNumber) {
        for (Map.Entry<String, Long> entry : refIndexMap.entrySet()) {
            if (chrNumber == entry.getValue()) {
                return entry.getKey();
            }
        }
        return null;
    }

    public String checkRefIndex() {
        
        if(this.refIndex.size()>0){
            
            return "refIndex Has value : EX key 0 is "+this.refIndex.get(1);
        }else{
            return "refIndex Has No value";    
        }
    }

    public String generateReadme(String svType){
        String readme = "";
        if(svType.equals("TD")){
            readme = "Tandem Report Format\n\n"
                    + ">Group\tGroup ID\tChromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tCorrectness\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back\n"
                    + "Extend view of junction [Reference]\n"
                    + "Support reads sequence [number of line equal to number of support reads]\n"
                    + ".\n"
                    + ".\n"
                    + ".\n"
                    + ">[Next Group]\n\n"
                    + "/********************************************/\n"
                    + "Junction can split into two part front and back\n"
                    + "Strand code : 0 = strand + and 1 = strand -\n"
                    + "For tandem Correctness higher is better";
        }else if(svType.equals("D")){
            readme = "Deletion Report Format\n\n"
                    + ">Group\tGroup ID\tChromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tCorrectness\tDeletion Size\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back\n"
                    + "Extend view of junction [Reference]\n"
                    + "Support reads sequence [number of line equal to number of support reads]\n"
                    + ".\n"
                    + ".\n"
                    + ".\n"
                    + ">[Next Group]\n\n"
                    + "/********************************************/\n"
                    + "Junction can split into two part front and back\n"
                    + "Strand code : 0 = strand + and 1 = strand -\n"
                    + "For deletion Correctness lower is better";
        }else if(svType.equals("IE")){
            readme = "InterInsertion Report Format\n\n"
                    + ">Group\tGroup ID\t[Front Junction] Chromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back"
                    + "\t||\t"
                    + "[Back Junction] Chromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back"
                    + "\tInsert Portion Size\tInsertion Junction\n"
                    + "Extend view of junction [Reference]\n"
                    + "Support reads sequence of both junction[number of line equal to number of support reads]\n"
                    + ".\n"
                    + ".\n"
                    + ".\n"
                    + ">[Next Group]\n\n"
                    + "/********************************************/\n"
                    + "Junction can split into two part front and back\n"
                    + "Strand code : 0 = strand + and 1 = strand -";                    
        }else if(svType.equals("IA")){
            readme = "IntraInsertion Report Format\n\n"
                    + ">Group\tGroup ID\t[Front Junction] Chromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back"
                    + "\t||\t"
                    + "[Back Junction] Chromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back"
                    + "\tInsert Portion Size\tInsertion Junction\n"
                    + "Extend view of junction [Reference]\n"
                    + "Support reads sequence of both junction[number of line equal to number of support reads]\n"
                    + ".\n"
                    + ".\n"
                    + ".\n"
                    + ">[Next Group]\n\n"
                    + "/********************************************/\n"
                    + "Junction can split into two part front and back\n"
                    + "Strand code : 0 = strand + and 1 = strand -";                    
        }else if(svType.equals("CH")){
            readme = "Chimeric Report Format\n\n"
                    + ">Group\tGroup ID\tChromosome Front:Position Front:Strand Front\tChromosome Back:Position Back:Strand Back\tSupport Reads\tNumber of Strand Pattern\tAnnotation Front\tAnnotation Back\n"
                    + "Extend view of junction [Reference]\n"
                    + "Support reads sequence [number of line equal to number of support reads]\n"
                    + ".\n"
                    + ".\n"
                    + ".\n"
                    + ">[Next Group]\n\n"
                    + "/********************************************/\n"
                    + "Junction can split into two part front and back\n"
                    + "Strand code : 0 = strand + and 1 = strand -\n"
                    + "For deletion Correctness lower is better";
        }
        
        return readme;
    }
    
    public void FilterIntraInsertion(){
        /**
         * this function will filter intrainsertionList by filter out SVPair that has a combination of DEl DEl out and save to separate list (intraInsertionList_outFilter)
         * It will replace the intrainsertionList with new set of filtered SVPair.
         */
        ArrayList<SVGroupPair> intraInsertionList_filter = new ArrayList(); 
        for(int i =0 ;i<this.intraInsertionList.size();i++){
            SVGroupPair svPair = this.intraInsertionList.get(i);
            if(svPair.istDelDel()){
                this.intraInsertionList_outFilter.add(svPair);
            }else{
                intraInsertionList_filter.add(svPair);
            }
        }
        
        this.intraInsertionList.clear();
        this.intraInsertionList.addAll(intraInsertionList_filter);
        
        intraInsertionList_filter.clear();
    }
}
