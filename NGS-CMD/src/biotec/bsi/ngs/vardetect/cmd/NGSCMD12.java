
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
import java.io.IOException;
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author worawich
 * 
 * Use for clustering test
 */
public class NGSCMD12 {

    public static void main(String[] args) throws IOException {
        
        String filename = "hg19_Format_AlignSortedCutResultMap(100k-200k)";
        AlignmentResultRead readAlign = SequenceUtil.readAlignmentReport("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/"+filename+".txt");
        
        String path = "/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/"+filename;
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
        System.out.println("********** Do Calculate Euclidient Distance *********");
        readAlign.calculateEuclidientdistance(); // Must have this order before clustering
//        VisualizeResult.visualizeDistanceTable(align);
//        System.out.println("********** Do Write Result DistanceTable *********");
//        readAlign.writeDistanceTableToPath(path, "txt");
        
        System.out.println("********** Do Clustering Group *********");
        System.out.println("Number of Read in consider : "+ readAlign.getResult().size());
        ArrayList<ClusterGroup> groupResult = Clustering.clusteringGroup(readAlign, 100);
        readAlign.addGroupReult(groupResult);
        System.out.println("********** Do Write Group Result *********");
        readAlign.writeClusterGroupToPath(path, "txt");
        
        System.out.println(" check number of group : " + groupResult.size());
        VisualizeResult.visualizeClusterGoup(groupResult);
        
        System.out.println("********** Do Reconstruct Sequence *********");
        readAlign.enableReconstruct();
        System.out.println("********** Do Write Pattern Report *********");
        readAlign.writePatternReport(path, "txt");
    
    }
}
