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
    
    public static void createReferenceFromContig(String[] inputPath) throws IOException{
        
        /**
         * This function will combine set of reference contig file into one single reference file read contig from fasta file
         * Recreate it as a reference like chromosome file in human 
         * The amount of base limit for each chromosome is 28bit 
         * Also create index to trace back to each scaffold 
         */
        
        int mask28Bit = 268435455 ;
        int sliding = 60;
        int window = 60;
        
        File inputFile = new File(inputPath[0]);
        String path = inputFile.getParent();
        
//        StringBuffer seq = new StringBuffer();
        Map<Long,String> index = new LinkedHashMap();
        ArrayList<String> listSeq = new ArrayList();
        ArrayList<Integer> chrList = new ArrayList();
        ArrayList<Integer> iniPosList = new ArrayList();
        ArrayList<Integer> lastPosList = new ArrayList();
        ArrayList<String> contigNameList = new ArrayList();
        ArrayList<String> refFileNameList = new ArrayList();
//        int iniPos = 0;
//        int lastPos = 0;
        int chrCount = 0;
//        boolean firstFlag = true;
//        String scaffoldName = null;
        
        File refFile = new File(path,"ContigRef.fa"); //File object
        FileWriter writer;
        if(refFile.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(refFile,false);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(refFile);
        }


        for(int numF=0;numF<inputPath.length;numF++){
            chrCount = chrCount+1;
            inputFile = new File(inputPath[numF]);     
            String name = inputFile.getName().split("\\.")[0];

            StringBuffer seq = new StringBuffer();
            int iniPos = 0;
            int lastPos = 0;
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
                            refFileNameList.add(name);

                            for(int i=0;i<1000;i++){
                               /**
                                * append "N" 1000 character
                                */
                               seq.append("N");
                            }

                            iniPos = seq.length();              // index is start at 0, the next iniPos is actual index plus 1 but to get actual index we have to minus 1. So, both cancel out and we don't have to do any thing.

                        }

                        scaffoldName = sb.substring(1);
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
                refFileNameList.add(name);
                listSeq.add(seq.toString());
                /*******************************/

            }catch(IOException e){

            }
        }
        
        
//        File refFile = new File(path,"ContigRef.fa"); //File object
//        FileWriter writer;
//        if(refFile.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(refFile,false);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(refFile);
//        }
        
        /**
         * Write reference file
         */
        
//        for(int i=0;i<listSeq.size();i++){
//            int numChr = i+1;
//            writer.write(">chr"+numChr+"\n");
//            String dummyString = listSeq.get(i);
//            writer.write(listSeq.get(i));           
//        }
//        writer.flush();
//        writer.close();
        for(int i=0;i<listSeq.size();i++){
            int numChr = i+1;
            writer.write(">chr"+numChr+"\n");
            String dummyString = listSeq.get(i);
            int n = (dummyString.length()-window)/sliding;       
            for(int numWin=0;numWin<n;numWin++){
                String s = dummyString.substring(numWin*sliding,numWin*sliding+window);
                writer.write(s+"\n");
            }
        }
        writer.flush();
        writer.close();
        
        
        
        
        
        
        File f = new File(path,"ContigRef.index"); //File object
        FileWriter writerIdx;
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writerIdx = new FileWriter(f,false);
        }else{
//            ps = new PrintStream(filename);
            writerIdx = new FileWriter(f);
        }

        /**
         * write index file
         */

        for(int i=0;i<chrList.size();i++){
            writerIdx.write(String.format("%s,%s,chr%d,%d,%d", refFileNameList.get(i),contigNameList.get(i),chrList.get(i),iniPosList.get(i),lastPosList.get(i)));
            writerIdx.write("\n");
        }
        writerIdx.flush();
        writerIdx.close();

    }

}