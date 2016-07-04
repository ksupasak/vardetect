/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author worawich
 */
public class AlignmentResultRead {
    private ArrayList<ShortgunSequence> shrtRead;
    private ArrayList<ClusterGroup> clusterResult;
    long[] allClusterCodeSorted;
    long[] allClusterCode;
    private static long mask = 268435455;
    
    public AlignmentResultRead(){
        
       this.shrtRead = new ArrayList(); 
    }
    
    public void addResult(ShortgunSequence inRead){
        this.shrtRead.add(inRead);
    }
    
    public ArrayList<ShortgunSequence> getResult(){
    
        return this.shrtRead;
    }
    
    public void createGroupingResult(){
        long dummyCode = 0;
        long oldDummyCode = 0;
        
        ClusterGroup group = new ClusterGroup();
        for(int i =0;i<this.allClusterCodeSorted.length;i++){
            dummyCode = this.allClusterCodeSorted[i];
            for(int j =0;j<this.shrtRead.size();j++){
                ShortgunSequence dummySS = shrtRead.get(j);
                if(dummySS.getClusterCode() == dummyCode){
                    if(i == 0){
                        group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }else if(i!=0 && Math.abs(dummyCode-oldDummyCode)<=100){
                        group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }else if(i!=0 && Math.abs(dummyCode-oldDummyCode)>100){
                        group = new ClusterGroup();
                        group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }
                    
                    System.out.println(dummyCode + "\t" + dummySS.getReadName());
                }
            }
        }
    }
    
    public void createAllClusterCode(){
        this.allClusterCode = new long[this.shrtRead.size()];
        for(int i =0;i<this.shrtRead.size();i++){
            long dummyCode = shrtRead.get(i).getClusterCode();
            this.allClusterCode[i] = dummyCode;
        }
    }
    
    public void createAllClusterCodeSorted(){
        this.allClusterCodeSorted = new long[this.shrtRead.size()];
        for(int i =0;i<this.shrtRead.size();i++){
            long dummyCode = shrtRead.get(i).getClusterCode();
            this.allClusterCodeSorted[i] = dummyCode;
        }
        Arrays.sort(this.allClusterCodeSorted);
    }
    
    public long[] getAllClusterCode(){
        
        //Arrays.sort(this.allClusterCode);
        return this.allClusterCode;
    }
    
    public long[] getAllClusterCodeSorted(){
        
        return this.allClusterCodeSorted;
    }
    
    public void writeSortedResultToPath(String path, String fa) throws FileNotFoundException, IOException {

       
        PrintStream ps = new PrintStream(path+"_AlignSortedResult."+ fa);
        
        
        for (int i=0;i<this.shrtRead.size();i++){           // Loop Mer by Mer
            
            ShortgunSequence dummySS = this.shrtRead.get(i);
        //--------------------------    

        //---------------------------------
            
            Map<Long,long[]> countMap =  dummySS.getAlignmentCountSorted();
            ps.println(">Alignment result of "+ dummySS.getReadName());
            ps.printf("%-30s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>28;
                long[] numCountPlusColor = countMap.get(positionCode);
                long numCount = numCountPlusColor[0];
                long red = numCountPlusColor[1];
                long yellow = numCountPlusColor[2];
                long orange = numCountPlusColor[3];
                long green = numCountPlusColor[4];
                long redInt = numCountPlusColor[5];
                long yellowInt = numCountPlusColor[6];
                long orangeInt = numCountPlusColor[7];
                long greenInt = numCountPlusColor[8];
                
                
                ps.format("Chr %d : Position %d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
            }
            ps.println();
        }
    }

    public void writeUnSortedResultToPath(String path, String fa) throws FileNotFoundException, IOException {

       
        PrintStream ps = new PrintStream(path+"_AlignUnSortedResult."+ fa);
        
        
        for (int i=0;i<this.shrtRead.size();i++){           // Loop Mer by Mer
            
            ShortgunSequence dummySS = this.shrtRead.get(i);
        //--------------------------    

        //---------------------------------
            
            Map<Long,long[]> countMap =  dummySS.getAlignmentCount();
            ps.println(">Alignment result of "+ dummySS.getReadName());
            ps.printf("%-30s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>28;
                long[] numCountPlusColor = countMap.get(positionCode);
                long numCount = numCountPlusColor[0];
                long red = numCountPlusColor[1];
                long yellow = numCountPlusColor[2];
                long orange = numCountPlusColor[3];
                long green = numCountPlusColor[4];
                long redInt = numCountPlusColor[5];
                long yellowInt = numCountPlusColor[6];
                long orangeInt = numCountPlusColor[7];
                long greenInt = numCountPlusColor[8];
                
                
                ps.format("Chr %d : Position %d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
            }
            ps.println();
        }
    }

    public void writeSortedCutResultToPath(String path, String fa, int threshold) throws FileNotFoundException, IOException {

        /* Must specify threshold for cut result (The result that less than threshold will be cut out)*/
        PrintStream ps = new PrintStream(path+"_AlignSortedCutResult."+ fa);
        
        
        for (int i=0;i<this.shrtRead.size();i++){           // Loop Mer by Mer
            
            ShortgunSequence dummySS = this.shrtRead.get(i);
        //--------------------------    

        //---------------------------------
            
            Map<Long,long[]> countMap =  dummySS.getAlignmentCountSortedCut(threshold);
            ps.println(">Alignment result of "+ dummySS.getReadName());
            ps.printf("%-30s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>28;
                long[] numCountPlusColor = countMap.get(positionCode);
                long numCount = numCountPlusColor[0];
                long red = numCountPlusColor[1];
                long yellow = numCountPlusColor[2];
                long orange = numCountPlusColor[3];
                long green = numCountPlusColor[4];
                long redInt = numCountPlusColor[5];
                long yellowInt = numCountPlusColor[6];
                long orangeInt = numCountPlusColor[7];
                long greenInt = numCountPlusColor[8];
                
                
                ps.format("Chr %d : Position %d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
            }
            ps.println();
        }
    }    
        
    
}
