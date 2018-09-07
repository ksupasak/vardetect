/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author worawich
 */
public class SumVCFresult {
    public static void main (String args[]) throws IOException{
        /**
         * create sum file of all sv type
         */
        Charset charset = Charset.forName("US-ASCII");
        Boolean firstFlag = true;
        String vcfType = "diploidSV";
        /************/
        
        /**
         * Modified map to new file name
         */
        Map<String,String> nameMap = new HashMap();
        String name_file = "/Volumes/4TB_WD/TB/Martin_new_list/TB_remove_dup_bam/Sample_fullName.txt";
        File nameFile = new File(name_file);
        Path pathF = Paths.get(nameFile.getPath());
        try (BufferedReader reader = Files.newBufferedReader(pathF, charset)) {
            String line = null;    
            while ((line = reader.readLine()) != null) {
                String[] portion = line.split("_");
                nameMap.put(portion[0], line);
            }
        }
        
        /**************/
        
        /**
         * start loop through dir
         */
        String myDirectoryPath = "/Volumes/4TB_WD/TB/Martin_new_list/TB_remove_dup_bam/manta_result/vcf";
        
        FileWriter writer;

        File dir = new File(myDirectoryPath);
        File[] directoryListing = dir.listFiles();
        if (directoryListing != null) {
            // loop in initial path dir
            for (File child : directoryListing) {
                // initial path folder level
                String[] fileComp = child.getName().split("\\.");
                if(fileComp[fileComp.length-1].equals("gz") && fileComp[fileComp.length-2].equals("vcf") && fileComp[fileComp.length-3].split("_")[1].equals(vcfType)){
                    String fileName = child.getName();
                    String sampleName = fileComp[fileComp.length-3].split("_")[0];
                    /**
                    * Modified map to new file name
                    */
                    sampleName = nameMap.get("ERR"+sampleName);
                    /**************/
                    File sumFile = new File(myDirectoryPath + File.separator + "sumResult_vcf.txt"); //File object        
                    if(sumFile.exists()){
//                                        ps = new PrintStream(new FileOutputStream(filename,true));
                        writer = new FileWriter(sumFile,true);
                        firstFlag = false;
//                                        writer = new FileWriter(summaryReport); //not append
                    }else{
//                                        ps = new PrintStream(filename);
                        writer = new FileWriter(sumFile);
                        firstFlag = true;
                    }

                    
                    InputStream fileStream = new FileInputStream(child);
                    InputStream gzipStream = new GZIPInputStream(fileStream);
                    Reader decoder = new InputStreamReader(gzipStream, charset);
//                    BufferedReader buffered = new BufferedReader(decoder);
                    
                    
                    Path path = Paths.get(child.getPath());
//                    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
                    try (BufferedReader reader = new BufferedReader(decoder)){
                        String line = null;
                        String data = null;
                        int count = 0;
                        
                        while ((line = reader.readLine()) != null) {
                            
                            if(line.charAt(0) == '#' && line.charAt(1)== '#'){
                                continue;
                            }else if(line.charAt(0) == '#' && line.charAt(1)!= '#'){
                                data = line;
                            }else if(line.charAt(0) != '#'){
                                data = line;
                            }
                            
                            
                            if(firstFlag == true && count == 0){

                                String modifLine = "SampleName\t"+data;
                                writer.write(modifLine);
                                writer.write("\n");

                            }else if(count > 0){
                                String modifLine = sampleName+"\t"+data;
                                writer.write(modifLine);
                                writer.write("\n");
                            }

                            count++;
                        }
                    }
                    writer.flush();
                    writer.close();
                }
            }
        }
    }
}
