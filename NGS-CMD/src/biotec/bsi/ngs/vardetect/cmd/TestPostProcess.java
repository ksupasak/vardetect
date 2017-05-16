/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.VariationResult;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import java.io.File;
import java.io.IOException;


/**
 *
 * @author worawich
 */
public class TestPostProcess {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/dm6_O4_4thread_th5_alignmentResult_longRead_part";
        String indexFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/drosophila/d.melanogaster/dm6_filter.index";            // use for traceback to natural chromosome name
        String saveFileType = "txt";
        int numPart = 1;
//        int readLength = 24;
        int merLength = 18;
        int maxFullMatch = 100; // it's percentage value
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/simLongRead/";
        
        
        File mainFile = new File(filename);
        String path = mainFile.getPath();
        String saveFileName = mainFile.getName().split("part")[0]+"_forLinuxSort";
//        File pathVar = new File(args[1]);
//        filename = pathVar.getName().split("\\.")[0];
//        path = pathVar.getParent()+"/";
//        saveFileName = pathVar.getPath().split("\\.")[0];
        
        String saveFileNameERR = "hg38_simData_mul_alignmentResult_th5_ERR_forLinuxSortV2";
        String saveFileName3661 = "hg38_simData_mul_alignmentResult_th5_forLinuxSort";
        String saveFileNameForPostProcess = "hg38_simData_mul_alignmentResult_th5_PostLinuxSortedV2";
        String filenameRNA = "hg38_Tha3.7_alignmentResult_part";
        String saveFilenameRNA = "hg38_Tha3.7_alignmentResult_forLinuxSort";
        String filenameERR = "hg38_err1_alignmentResult_part";
        String saveFilenameERR = "hg38_err1_alignmentResult_forLinuxSort_part";
        
        
        
        for(int i=1;i<=numPart;i++){
            
            AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2(filename+i+".txt",merLength);
//            AlignmentResultRead readAlign = SequenceUtil.readBinaryAlignmentReportV3(filename+i+".bin",merLength);
            
            readAlign.countAlignMatch();
            readAlign.writeMatchCountReport(filename+i,indexFile);              // write bed graph file format for overall vistualize
            
            System.out.println("Begin create color array");
            Clustering.createColorArrayV2(readAlign, merLength);        
            System.out.println("Done create color array");
            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType,"all",maxFullMatch);
//            VariationResult varRes = SequenceUtil.analysisResultFromFile(path+saveFileNameForPostProcess+".txt",18,100);
            readAlign = null;
            System.gc();
        }
    }
  
}
