/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestTranslateData {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/simLongRead/hg38_simLong_4thread_th5_alignmentResult_LongRead_part";
        String indexFile = "";
        String inputDataType = "bin";
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
        
        AlignmentResultRead readAlign = null;
        
        for(int i=1;i<=numPart;i++){
            
            if(inputDataType == "bin"){
                readAlign = SequenceUtil.readBinaryAlignmentReportV3(filename+i+".bin",merLength);
            }else if(inputDataType == "txt"){
                readAlign = SequenceUtil.readAlignmentReportV2(path+filename+i+".txt",merLength);
            }

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
