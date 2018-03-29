/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 *
 * @author worawich
 */
public class testExecuteCommandLine {
    
    public static void main(String[] args) throws IOException{
        
        String bamFile = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/rmDup_Bam/NA12878-L1.rmdup.bam";
        String chrName = "8";
        String breakPoint = "58123049";
        
        
        String baseCommand = "samtools view " + bamFile;
        
        
        String addCommand = "chr"+chrName+":"+breakPoint+"-"+breakPoint;
        String command = baseCommand + " " + addCommand;
//        String command = baseCommand + " " + addCommand + " | wc";
//        String domainName = "google.com";
		
        //in mac oxs
//        command = "ping -c 3 " + domainName;


        command = "/usr/local/bin/samtools view /Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/rmDup_Bam/NA12878-L1.rmdup.bam chr8:58123049-58123049 | wc";
        String resultF = executeCommandLine(command);
        String[] splitResF = resultF.split("\\s+");
        int coverageF = Integer.parseInt(splitResF[1]);
        
    }
    
    public static String executeCommandLine(String command) throws IOException {

            StringBuffer output = new StringBuffer();
            String[] cmd = {"/bin/sh","-c",command};
            
//            ProcessBuilder pb = new ProcessBuilder("/bin/sh", "-c", command);
//            Process p = pb.start();
//            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
//            String line = "";			
//            while ((line = reader.readLine())!= null) {
//                    output.append(line + "\n");
//            }
            Process p;

            try {
                    p = Runtime.getRuntime().exec(cmd);
//                    p = Runtime.getRuntime().exec(command);
                    p.waitFor();
                    BufferedReader reader = 
                        new BufferedReader(new InputStreamReader(p.getInputStream()));

                    String line = "";			
                    while ((line = reader.readLine())!= null) {
                            output.append(line + "\n");
                    }
//                    p.waitFor();
            } catch (Exception e) {
                    e.printStackTrace();
            }

            return output.toString();

    }
    
//    Process linuxControlProcess;
//        try{
//            String[] cmd = {"/bin/sh","-c","sort -t, -k10,10 -k9,9n "+fullPathSaveUnSortFile+" >> "+fullPathSaveSortFile};
//            linuxControlProcess = Runtime.getRuntime().exec(cmd);
//            linuxControlProcess.waitFor();
//            String[] cmd2 = {"/bin/sh","-c","rm "+fullPathSaveUnSortFile};
//            linuxControlProcess = Runtime.getRuntime().exec(cmd2);
//            linuxControlProcess.waitFor();
//        } catch (Exception e){
//            e.printStackTrace();
//        }
    
}
