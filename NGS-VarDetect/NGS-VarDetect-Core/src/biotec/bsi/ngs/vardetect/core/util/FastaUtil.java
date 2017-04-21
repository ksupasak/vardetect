/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;



/**
 *
 * @author worawich
 */
public class FastaUtil {
    
    public static void createReferenceFromContig(String inputPath) throws IOException{
        
        /**
         * This function will read contig from fasta file
         * Recreate it as a reference like chromosome file in human 
         * The amount of base limit for each chromosome is 28bit 
         * Also create index to trace back to each scaffold 
         */
        
        int mask28Bit = 268435455 ;
        
        
        File inputFile = new File(inputPath);     
        
        
        
        StringBuffer seq = new StringBuffer();
        Map<Long,String> index = new LinkedHashMap();
        ArrayList<String> listSeq = new ArrayList();
        ArrayList<Integer> chrList = new ArrayList();
        ArrayList<Integer> iniPosList = new ArrayList();
        ArrayList<Integer> lastPosList = new ArrayList();
        ArrayList<String> contigNameList = new ArrayList();
        int iniPos = 0;
        int lastPos = 0;
        int chrCount = 1;
        boolean firstFlag = true;
        String scaffoldName = null;
        
        
        /**
         * Read file
         */
        try(BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(inputFile)))){
            String sb = null;
            while ((sb = reader.readLine()) != null){
                
                if(sb.charAt(0)=='>'){
                    
                    
                    if(firstFlag == false){
                        lastPos = seq.length()-1;           // index is start at 0 So, we have to minus 1 in order to get actual index
                        
//                        long indexCode = ((long)chrCount<<28|(long)iniPos)<<28|(long)lastPos;       // iniPos and lastPos is actual index. it can be used directly
//                        
//                        index.put(indexCode,scaffoldName);
                        
                        chrList.add(chrCount);
                        iniPosList.add(iniPos);
                        lastPosList.add(lastPos);
                        contigNameList.add(scaffoldName);
                        
                        for(int i=0;i<1000;i++){
                           /**
                            * append "N" 1000 character
                            */
                           seq.append("N");
                        }
                        
                        iniPos = seq.length();              // index is start at 0, the next iniPos is actual index plus 1 but to get actual index we have to minus 1. So, both cancel out and we don't have to do any thing.
                        
                    }
                    
                    scaffoldName = sb;
                    firstFlag = false;
                    
                }else{
                    
                    if(sb.length()>18){
                        
                        if(seq.length()+sb.trim().length() > mask28Bit){
                            
                            /**
                             * Reach maximum number of base in one seq
                             * initiate new seq and also reset iniPos to zero and firstFlag to true
                             * Reset or Update: 
                             *  iniPos to zero 
                             *  firstFlag to true
                             *  chrCount plus one
                             */
                            
                            listSeq.add(seq.toString());
                            seq = new StringBuffer();
                            seq.append(sb.trim());
                            chrCount++;
//                            firstFlag = true;
                            iniPos = 0;                            
                        }else{
                            seq.append(sb.trim());
                        }                        
                    }else{
                        firstFlag = true;
                    }
                } 
            }
            
            /**
             * for last set
             */
            lastPos = seq.length()-1;           // index is start at 0 So, we have to minus 1 in order to get actual index

            chrList.add(chrCount);
            iniPosList.add(iniPos);
            lastPosList.add(lastPos);
            contigNameList.add(scaffoldName);
            listSeq.add(seq.toString());
            /*******************************/
            
        }catch(IOException e){
            
        }
        
        /**
         * Write file
         */
        
        for(int i=0;i<listSeq.size();i++){
            int numFile = i+1;
            String[] dummyS = inputPath.split(".fasta");
            String fileName = dummyS[0]+"_"+numFile+".fa";
            File f = new File(fileName); //File object
            FileWriter writer;
            if(f.exists()){
    //            ps = new PrintStream(new FileOutputStream(filename,true));
                writer = new FileWriter(f,true);
            }else{
    //            ps = new PrintStream(filename);
                writer = new FileWriter(f);
            }

            writer.write(listSeq.get(numFile));
            writer.flush();
            writer.close();
        }
        
        
        /**
         * write index file
         */
        String fileName = inputPath.split(".")[0]+".index";
        File f = new File(fileName); //File object
        FileWriter writer;
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(f,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(f);
        }
            
        for(int i=0;i<chrList.size();i++){
            writer.write(String.format("%s,%d,%d,%d", contigNameList.get(i),chrList.get(i),iniPosList.get(i),lastPosList.get(i)));
            writer.write("\n");
        }
        writer.flush();
        writer.close();
        
    }

}
