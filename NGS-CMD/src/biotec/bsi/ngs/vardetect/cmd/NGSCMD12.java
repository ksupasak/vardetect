
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author worawich
 * 
 * Use for clustering test (local map implementation)
 */
public class NGSCMD12 {

    public static void main(String[] args) throws IOException {
        
//        String filename = "hg19_Format_AlignSortedCutResultMap(100k-200k)";
        String inputFileName = "/Users/worawich/VMdev/Siriraj/JT/JT.unmapped.sam";
        String filename = "hg19JT_unalign_Format_AlignSortedCutResultMap_part1";
        AlignmentResultRead readAlign = SequenceUtil.readAlignmentReport("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/"+filename+".txt");
        ArrayList<AlignmentResultRead> localAlignRes = new ArrayList<AlignmentResultRead>();
        Map<Integer,ArrayList<String>> groupMap = new HashMap<Integer,ArrayList<String>>();
        ArrayList<Map<Integer,ArrayList<String>>> groupList = new ArrayList<Map<Integer,ArrayList<String>>>();
        
        String path = "/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/"+filename;
        String savePath = "/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/";
        String saveFilename = filename+"_clusterGroup_th74";
        /* Old Grouping algorithm
        align.createAllClusterCode();
        align.createAllClusterCodeSorted();
        
        for(int i =0;i<align.getAllClusterCode().length;i++){
            System.out.println("Check cluster code: " + align.getAllClusterCode()[i]);
        }
        for(int i =0;i<align.getAllClusterCode().length;i++){
            System.out.println("Check cluster code sorted: " + align.getAllClusterCodeSorted()[i]);
        }
        
        align.createGroupingResult();
        System.out.println(" ****** Check cluster result: " + align.getclusterResult().size());
        */
        Map<Long,ArrayList<String>> chrMatchMap = Clustering.createChrMatchMap(readAlign);
        Set set = chrMatchMap.keySet();
        Iterator chrNumberIter = set.iterator();
        
        System.out.println("Try to class Group");
        int num = 0;
        while(chrNumberIter.hasNext()){                                                     // Loop over chromosome group which has 24 group
            ArrayList<String> readNameList = chrMatchMap.get(chrNumberIter.next());
            InputSequence inSeq = SequenceUtil.readSamFile(inputFileName, readNameList);
            groupMap = SequenceUtil.localAlignment(inSeq, 18, 1, 74);                      // inside localAlignment will have loop over member in this chromosome group (all possible pair of each member)
            if(groupMap.isEmpty()){
                //groupList.add(0);
            }else{
                groupList.add(groupMap);
            }
            num++;
            System.out.println("Done "+num+"/24");
        } 
        /* 
        Do local alignment between readin the same group in Map
        */
        System.out.println("Filter Group");
        int minCoverage = 2;
        Map<Integer,ArrayList<String>> groupResult = new HashMap<Integer,ArrayList<String>>();
        groupResult = Clustering.filterClusterGroupLocalAlignment(groupList, minCoverage);
        System.out.println("Save cluster Result");
        Clustering.writeLocalAlignmentInFile(groupResult, savePath, saveFilename);
        System.out.println("Done");
//        System.out.println("********** Do Calculate Euclidient Distance *********");
//        readAlign.calculateEuclidientdistance(); // Must have this order before clustering
////        VisualizeResult.visualizeDistanceTable(align);
////        System.out.println("********** Do Write Result DistanceTable *********");
////        readAlign.writeDistanceTableToPath(path, "txt");
//        
//        System.out.println("********** Do Clustering Group *********");
//        System.out.println("Number of Read in consider : "+ readAlign.getResult().size());
//        ArrayList<ClusterGroup> groupResult = Clustering.clusteringGroup(readAlign, 100);
//        readAlign.addGroupReult(groupResult);
//        System.out.println("********** Do Write Group Result *********");
//        readAlign.writeClusterGroupToPath(path, "txt");
//        
//        System.out.println(" check number of group : " + groupResult.size());
//        VisualizeResult.visualizeClusterGoup(groupResult);
//        
//        System.out.println("********** Do Reconstruct Sequence *********");
//        readAlign.enableReconstruct();
//        System.out.println("********** Do Write Pattern Report *********");
//        readAlign.writePatternReport(path, "txt");
    
    }
}
