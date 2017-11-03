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
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author worawich
 */
public class CreateAlignedSample {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Error/ERR718196_1_alnRes_part";
        String indexFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/drosophila/d.melanogaster/dm6_filter.index";            // use for traceback to natural chromosome name
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/O3_S3_L001_R2_001.fa";
        String saveFileType = "txt";
        int numPart = 4;
//        int readLength = 24;
        int merLength = 12;
        int maxFullMatch = 100; // it's percentage value
        int propotion = 1000000;
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/simLongRead/";
        long targetFront = 4092077;
        long targetBack = 4092919;
        int targetBoundary = 500;
        
        File mainFile = new File(filename);
        String path = mainFile.getParent();
        String saveFileName = mainFile.getName().split("part")[0]+"forLinuxSortTarget";
        String saveBedFileName = mainFile.getName().split("part")[0];
        String saveSampleFileName = mainFile.getName().split("part")[0]+"Sample";
//        File pathVar = new File(args[1]);
//        filename = pathVar.getName().split("\\.")[0];
//        path = pathVar.getParent()+"/";
//        saveFileName = pathVar.getPath().split("\\.")[0];

        AlignmentResultRead readAlign = new AlignmentResultRead();   
        for(int i=1;i<=numPart;i++){
            
            readAlign = SequenceUtil.readAlignmentReportV2(filename+i+".txt",merLength);
//            readAlign = SequenceUtil.readAlignmentReportFilterToTarget(readAlign,filename+i+".txt", merLength, targetFront, targetBack, targetBoundary);
//            AlignmentResultRead readAlign = SequenceUtil.readBinaryAlignmentReportV3(filename+i+".bin",merLength);
//            readAlign.writeSortedCutResultMapToPathInFormatLongRead(folder.toString()+File.separator,savefilename, filetype);
            
            
//            System.out.println("Begin create color array");
//            Clustering.createColorArrayV2(readAlign, merLength);        
//            System.out.println("Done create color array");
////            readAlign.sortFilterResult("all",maxFullMatch);
//            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType,"all",maxFullMatch);           
//            VariationResult varRes = SequenceUtil.analysisResultFromFile(path+saveFileNameForPostProcess+".txt",18,100);    
        }
        /**
        * Generate Bed graph  
        */
        readAlign.countAlignMatch();
        readAlign.writeMatchCountReport(filename);              // write bed graph file format for overall vistualize
//        readAlign.writeMatchCountReport(filename+i,indexFile);              // write bed graph file format for overall vistualize + index file
        /***************************************/
        
        System.out.println("Begin create color array");
        Clustering.createColorArrayV2(readAlign, merLength);        
        System.out.println("Done create color array");
//            readAlign.sortFilterResult("all",maxFullMatch);
        readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType,"all",maxFullMatch);
        readAlign.writeAlignmentSummaryReport(filename);
//        *****for create set of sample that align********//            
//            readAlign.writeAlignSequenceReadFasta(path, saveSampleFileName, sampleFile, propotion, "all",maxFullMatch);
//*****************//

        readAlign = null;
        System.gc();
    }
  
}
