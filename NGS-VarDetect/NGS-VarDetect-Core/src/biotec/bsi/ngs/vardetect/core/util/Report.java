/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class Report {
    
    public static void readCoverageReport(String fullpath){
        File file = new File(fullpath);
        String[] filenameComponent = file.getName().split("\\.")[0].split("_");
        int numFilenameComponent = filenameComponent.length;
        String fileType = file.getName().split("\\.")[1];
        String variantFileType = filenameComponent[numFilenameComponent-1];
        String reportType = filenameComponent[numFilenameComponent-2];
        
    }
    
    public static void reOrderGroupNumber(String filename) throws IOException{
        String[] fileNameComponenet = filename.split("\\.");
        String savefilename = fileNameComponenet[0]+"_reGroup.txt";
        FileWriter writer;        
        /**
         * Check File existing
         */
        
        File f = new File(savefilename); //File object        
        if(f.exists()){
            // append to file if file exist
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(savefilename,true);
        }else{
            // create new file if file not exist
//            ps = new PrintStream(filename);
            writer = new FileWriter(savefilename);
        }

        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(filename);
//        String name = null;
//        int actStart = readStart*2;     //this is actual start of line in file (compatible only specific file 3661 and 3662 .fasta file)
//        int actStop = readLimit*2;
    //    String seq = "";
        
        StringBuffer seq = new StringBuffer();
        
        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;
            String[] data = null;
            int groupCount = 0;
            while ((line = reader.readLine()) != null) {
                String newLine = "";
                if(line.charAt(0)=='G'){
                    groupCount++;
                    String groupHeadLine = line.toString();
                    String[] groupHeadLineComp = groupHeadLine.split("\t");
                    
                    groupHeadLineComp[0]="Group "+Integer.toString(groupCount);
                    
                    for(int i=0;i<groupHeadLineComp.length;i++){
                        newLine = newLine+groupHeadLineComp[i]+"\t";
                    }
//                    if(groupHeadLine.charAt(7) == ' '){
//                        newLine = 
//                        newLine = groupHeadLine.substring(0,6)+groupCount+groupHeadLine.substring(7);
//                    }else{
//                        newLine = groupHeadLine.substring(0,6)+groupCount+groupHeadLine.substring(8);
//                    }
                    
                    writer.write(newLine);
                    writer.write("\n");
                     
                }else{
                    writer.write(line);
                    writer.write("\n");              
                } 
            }
        }
        writer.flush();
        writer.close();
    }
    
}
