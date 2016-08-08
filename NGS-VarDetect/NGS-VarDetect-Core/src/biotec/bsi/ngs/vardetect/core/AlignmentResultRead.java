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
    private long[] allClusterCodeSorted;
    private long[] allClusterCode;
    private double[][] distanceTable;
    private static long mask = 268435455;
    private ClusterGroup group;
    
    public AlignmentResultRead(){
        
       this.shrtRead = new ArrayList();
       this.group = new ClusterGroup();
       this.clusterResult = new ArrayList();
    }
    
    public void addResult(ShortgunSequence inRead){
        this.shrtRead.add(inRead);
    }
    
    public void addGroupReult(ArrayList<ClusterGroup> inGroupResult){
        this.clusterResult = inGroupResult;
    }
    
    public ArrayList<ShortgunSequence> getResult(){
    
        return this.shrtRead;
    }
    
//    public void createdistancetable(){
//        int sizeMax = 2;
//        for(int i=0;i<this.shrtRead.size();i++){
//            ShortgunSequence dummyMainSS = shrtRead.get(i);
//            int sizeRes = dummyMainSS.countResultSortedCut.size();
//            if (sizeRes < sizeMax){
//                for (int k=0;k<sizeRes;k++){
//                    long dummyPos = (long)dummyMainSS.getListPosMatch().get(k);
//                    
//                }
//            }else if (sizeRes >= sizeMax){
//                for (int k=0;k<sizeRes;k++){
//                    long dummyPos = (long)dummyMainSS.getListPosMatch().get(k);
//                }
//            }
//            
//            
//            
//            for (int j=0;j<this.shrtRead.size();j++){
//                ShortgunSequence dummySubSS = shrtRead.get(j);
//                
//                if(dummyMainSS.getReadName() != dummySubSS.getReadName()){
//                    ArrayList dummyCheckSS = dummyMainSS.getListChrMatch();
//                    dummyCheckSS.retainAll(shrtRead);
//                    dummyMainSS.getListChrMatch().indexOf(dummyCheckSS);
//                    
//                    //dummyMainSS.getReadName()
//                    //dummyMainSS.getListChrMatch().reta
//                            
//                }
//            }
//        }
//        
//    }
    
    public double[][] getDistanceTable(){
        return this.distanceTable;
    }
    
    public void enableReconstruct(){
        for(int i =0;i<shrtRead.size();i++){
            shrtRead.get(i).detectStrandPattern();
        }
    }
    
    public void createGroupCharacteristic(double threshold){
        /* Create grouping characteristic of each shortgun sequence(store in shrtRead) */
        for(int i=0;i<this.shrtRead.size();i++){
            shrtRead.get(i).createInGroupOutGroup(threshold);
        }
    }
    
    public void calculateEuclidientdistance(){
        
        this.distanceTable = new double[shrtRead.size()][shrtRead.size()]; // Prelocate 2D double array for store distance value (size = [number of shortgun read] x [number of shortgun read])
        for(int i =0;i<shrtRead.size();i++){
            double[] distanceVector = new double[shrtRead.size()]; // distanceVector is a 1D double array stand for store distance value of main shortgun sequence clustervector and all other shortgun sequence clustervector
            ShortgunSequence dummyMainSS = shrtRead.get(i);
            for(int j=0;j<shrtRead.size();j++){
                ShortgunSequence dummySubSS = shrtRead.get(j);
                
                /*  distance function will calculate the distance vector of two in put vector               */
                /*  We will pass the clusterVector as a input for distance function                         */
                /*  Method getClusterVector will create the vector of align position of shortgun sequence   */
                /*  Ex: read0SS0 has align on chr1 pos: 1000 and chr21 pos: 200                             */
                /*  the clusterVector will look like this                                                   */
                /*  [1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,0,0,0]                                    */
                /*  as you can see the vector has 24 element which has value at first and 21 element follow by the chr that align */
                /*  In this Implementtation we selected only top two alignment result to create cluster vector */ 
                
                distanceVector[j] = distance(dummyMainSS.getClusterVector(),dummySubSS.getClusterVector()); // this will return single value of distance between two vector 
                this.distanceTable[i][j] = distanceVector[j];
//                System.out.println();
//                System.out.println("DummyMainSS/"+dummyMainSS.getReadName()+" pair with DummySubSS/"+dummySubSS.getReadName()+" : distanceVector["+j+"] = "+distanceVector[j]);
//                System.out.println();
            }
//            System.out.println("Check before add to "+dummyMainSS.getReadName());
            for(int a=0;a<distanceVector.length;a++){
//                System.out.print("\t"+distanceVector[a]);
            }
//            System.out.println();
            dummyMainSS.addDistanceVector(distanceVector); // store distance value in to main shrtgun sequence
//            System.out.println("Check after add to "+dummyMainSS.getReadName());
            for(int a=0;a<distanceVector.length;a++){
//                System.out.print("\t"+distanceVector[a]);
            }
//            System.out.println();
            
        }
//        for(int i =0;i<shrtRead.size();i++){
//            ShortgunSequence dummyMainSS = shrtRead.get(i);
//            System.out.println("************ ReadName:"+dummyMainSS.getReadName()+" check saved vector distance ************");
//            for (int check =0;check<shrtRead.size();check++){
//                System.out.print("\t"+dummyMainSS.getDistanceVector()[check]);   
//            }
//            System.out.println();
//        }
    }
    
    public double distance(long[] a, long[] b){
        double diff_square_sum = 0.0;
        for (int i = 0; i<a.length; i++){
            diff_square_sum += (a[i]-b[i]) * (a[i]-b[i]);
        }
        return Math.sqrt(diff_square_sum);
    }
    
    public void createGroupingResult(){
        long dummyCode = 0;
        long oldDummyCode = 0;
        
        this.group = new ClusterGroup();
        for(int i = 0;i<this.allClusterCodeSorted.length;i++){
            dummyCode = this.allClusterCodeSorted[i];
            for(int j =0;j<this.shrtRead.size();j++){
                ShortgunSequence dummySS = shrtRead.get(j);
                if(dummySS.getClusterCode() == dummyCode){
                    if(i == 0){
                        System.out.println(" Check : Do adding in first group (First time) i = " + i+ " : j =  "+j );
                        this.group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }else if(i!=0 && Math.abs(dummyCode-oldDummyCode)<=100){
                        System.out.println(" Check : Do adding in group : i = " + i+ " : j =  "+j);
                        this.group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }else if(i!=0 && Math.abs(dummyCode-oldDummyCode)>100){
                        this.clusterResult.add(this.group); // adding to array before renew it
                        System.out.println(" Check : Do create new group and add to new group : i = " + i+ " : j =  "+j);
                        
                        this.group = new ClusterGroup();
                        this.group.addShortgunRead(dummySS);
                        oldDummyCode = dummyCode;
                    }
                    
                    System.out.println(dummyCode + "\t" + dummySS.getReadName());
                }
            }
        }
        this.clusterResult.add(this.group); // adding to array (for last group)
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
    
    public ArrayList<ClusterGroup> getclusterResult(){
        
        return this.clusterResult;
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
            ps.printf("%-35s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","Strand","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>29; // 29 bit because positionCode have this structure [chr|strand|position]
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
                
                String strandNot = "no";
                if(((positionCode>>28)&1) == 1){
                    strandNot = "+";
                }else if(((positionCode>>28)&1) == 0){
                    strandNot = "-";
                }   
                
                
                ps.format("Chr %2d : Position %12d|\t%8s\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,strandNot,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
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
            ps.printf("%-35s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","Strand","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>29;
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
                
                String strandNot = "no";
                if(((positionCode>>28)&1) == 1){
                    strandNot = "+";
                }else if(((positionCode>>28)&1) == 0){
                    strandNot = "-";
                }   
                
                ps.format("Chr %2d : Position %12d|\t%8s\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,strandNot,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
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
//            ps.printf("%-30s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            ps.printf("%-35s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s%n","Result","Strand","NumMatch","Green","Yellow","Orange","Red","GreenInt","YellowInt","OrangeInt","RedInt");
            Set allPos = countMap.keySet();
            Iterator iterPos = allPos.iterator();
            while(iterPos.hasNext()){
                long positionCode = (long)iterPos.next();
                long alignPos = positionCode&mask;
                long chrNumber = positionCode>>29;
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
                
                String strandNot = "no";
                if(((positionCode>>28)&1) == 1){
                    strandNot = "+";
                }else if(((positionCode>>28)&1) == 0){
                    strandNot = "-";
                }   
                
//                ps.format("Chr %d : Position %d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
                ps.format("Chr %2d : Position %12d|\t%8s\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d%n",chrNumber,alignPos,strandNot,numCount,green,yellow,orange,red,greenInt,yellowInt,orangeInt,redInt);
            }
            ps.println();
        }
    }

    public void writeDistanceTableToPath(String path, String fa) throws FileNotFoundException, IOException {

       /* Must specify threshold for cut result (The result that less than threshold will be cut out)*/
        PrintStream ps = new PrintStream(path+"_DistanceTable."+ fa);
        
        ps.println("Distance Table");
        ps.format("Reads Name");
        for (int i=0;i<this.shrtRead.size();i++){           // Loop Mer by Mer
            
            ShortgunSequence dummySS = this.shrtRead.get(i);
        //--------------------------    

        //---------------------------------
 
            ps.format("\t%10s",dummySS.getReadName());
        }
        ps.println();
        for (int i=0;i<this.shrtRead.size();i++){ 
            ShortgunSequence dummySS = this.shrtRead.get(i);
            ps.print("Name: "+dummySS.getReadName());
            
            for(int j=0;j<this.shrtRead.size();j++){
                ps.format("\t%10.2f", dummySS.getDistanceVector()[j]);
            }
            ps.println();
        }   
    }    
        
    public void writeClusterGroupToPath(String path, String fa) throws FileNotFoundException, IOException {

       
        PrintStream ps = new PrintStream(path+"_ClusterGroup."+ fa);
       
        for(int i=0;i<this.clusterResult.size();i++){
            ClusterGroup dummyCG = this.clusterResult.get(i);
            ps.print("Group" + i + ":");
            ArrayList<ShortgunSequence> listSS = dummyCG.getShortgunRead();
            for(int j=0;j<listSS.size();j++){
                ShortgunSequence dummySS = listSS.get(j);
                ps.print("\t"+dummySS.getReadName());
            }
            ps.println();
        }  
    }
    
    public void writePatternReport(String path, String fa) throws FileNotFoundException, IOException {
        
        PrintStream ps = new PrintStream(path+"_PossiblePattern."+ fa);
        ps.println("Possible Fusion Pattern");
        ps.println();
        for(int i =0;i<this.shrtRead.size();i++){
        
            ShortgunSequence dummySS = this.shrtRead.get(i);
            ArrayList<ReconstructSequence> dummyListRecon = dummySS.getListReconSeq();
            ps.println("Read Name: " + dummySS.getReadName());
            
            for(int j=0;j<dummyListRecon.size();j++){
                ReconstructSequence dummyRecon = dummyListRecon.get(j);
                ps.println(""+dummyRecon.getResultString());
            }
            ps.println();
        }
        
    }
  
}
