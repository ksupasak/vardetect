/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author worawich
 */
public class NGSCMD_V2 {

    private static int minPeakF=10;
    private static int minPeakB=5;
    private static int merCov=50;
    private static int maxDup=1;
    private static int countDup=1;
    private static int numThread=1;
    private static long numSkip=0;
    private static int numRepeat=50;
    private static String refPath="";
    private static int numRead=10000;
    private static int numMer=16;
    private static int filterMode=3;
    
    public static void main(String[] args) throws IOException, FileNotFoundException, InterruptedException, ParseException {
        
        
        

        
        // create the command line parser
        CommandLineParser parser = new DefaultParser();

        // create the Options
        Options options = new Options();
//        options.addOption( "a", "all", false, "do not hide entries starting with ." );
//        options.addOption( "A", "almost-all", false, "do not list implied . and .." );
//        options.addOption( "b", "escape", false, "print octal escapes for nongraphic "
//                                                 + "characters" );
        options.addOption( Option.builder("r").longOpt("ref-file")
                                        .desc("Reference file path")
                                        .required()
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("t").longOpt("thread-num")
                                        .desc("number of Thread" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("n").longOpt("read-num")
                                        .desc("total number of read" )
                                        .required()
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("s").longOpt("skip-read")
                                        .desc("number of read to skip (use when debug)" )
                                        .hasArg()
                                        .argName("long")
                                        .build());
        options.addOption( Option.builder("d").longOpt("max-dup")
                                        .desc("maximum duplicate pattern" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("p").longOpt("mix-pattern")
                                        .desc("minimum peak pattern" )
                                        .hasArgs().numberOfArgs(2)
                                        .argName("int int")
                                        .build());
        options.addOption( Option.builder("m").longOpt("mer-size")
                                        .desc("size of mer" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("C").longOpt("mer-coverage")
                                        .desc("Filter value : number of mer coverage" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("F").longOpt("filter-mode")
                                        .desc("Filter mode has 1,2,3  number of mode" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("R").longOpt("ref-repeat")
                                        .desc("Filter value : number of base repeat" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("D").longOpt("filter-dup")
                                        .desc("Filter duplication count" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("h").longOpt("help")
                                        .desc("algorithm usage help" )
                                        .build());
        
//        options.addOption( "B", "ignore-backups", false, "do not list implied entried "
//                                                         + "ending with ~");
//        options.addOption( "c", false, "with -lt: sort by, and show, ctime (time of last " 
//                                       + "modification of file status information) with "
//                                       + "-l:show ctime and sort by name otherwise: sort "
//                                       + "by ctime" );
//        options.addOption( "C", false, "list entries by columns" );

//        args = new String[]{ "--block-size=10" };
//        args = new String[]{ "-A","-S 10" };

        try {
            // parse the command line arguments
            CommandLine line = parser.parse( options, args );

            
            if( line.hasOption( "p" ) ) {                
                System.out.println( line.getOptionValues("p")[0]);
                System.out.println( line.getOptionValues("p")[1]);
                minPeakF = Integer.parseInt(line.getOptionValues("p")[0]);
                minPeakB = Integer.parseInt(line.getOptionValues("p")[1]);
            }
            
            if(line.hasOption("C")) merCov = Integer.parseInt(line.getOptionValue("C"));
            if(line.hasOption("d")) maxDup = Integer.parseInt(line.getOptionValue("d"));
            if(line.hasOption("D")) countDup = Integer.parseInt(line.getOptionValue("D"));
            if(line.hasOption("F")) filterMode = Integer.parseInt(line.getOptionValue("F"));
            if(line.hasOption("m")) numMer = Integer.parseInt(line.getOptionValue("m"));
            if(line.hasOption("n")) numRead = Integer.parseInt(line.getOptionValue("n"));
            if(line.hasOption("r")) refPath = line.getOptionValue("r");
            if(line.hasOption("R")) numRepeat = Integer.parseInt(line.getOptionValue("R"));
            if(line.hasOption("s")) numSkip = Long.parseLong(line.getOptionValue("s"));
            if(line.hasOption("t")) numThread = Integer.parseInt(line.getOptionValue("t"));

            if(line.hasOption("h")){
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp( "SVP parameter setting help", options );
            }
        }
        catch( ParseException exp ) {
            System.out.println( "Unexpected exception:" + exp.getMessage() );
        }
        
        
        
        
        
        
//         *String refPath = args[0];
//         *CombineReferenceSequence ref = SequenceUtil.getCombineReferenceSequence(refPath,16); //runFile hg19.fa
//         
//         *ref.setMinimumPeakPattern(10, 5);

// for large batch with multi thread         
//         ref.setNumberOfThread(8);
//         ref.setTotalRead(86000000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(3);

// for large batch with multi thread         
//         ref.setNumberOfThread(4);
//         ref.setTotalRead(100000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(3);
         
// for small batch with multi thread         
//         *ref.setNumberOfThread(4);
//         ref.setTotalRead(10000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(1);
//         ref.setRandomAccess(true);
         

// for debug         
//         *ref.setNumberOfThread(1);
//         ref.setTotalRead(1);
//         ref.setSkipRead(3997);
//         ref.setMaximumDuplicatePattern(1);
//         ref.setMinimumPeakPattern(10, 2);
//         ref.setRandomAccess(true);
         
         
//         ref.searchMer(0);
//         ref.setOutputFile("/Users/soup/human/RB_cancer/277T_sorted_unmap.all.out");
//         ref.runProfileSV("/Users/soup/human/RB_cancer/277T_sorted_unmap.bam");
         
//         
//         
//         int found = 0;
//         for(int i=0;i<1000000;i++){
//         int v = (int)(Math.random()*Integer.MAX_VALUE);
//         int pos[] = ref.searchMer(v);
//         if(pos!=null){
//             found ++;
//            System.out.println(""+(i+1)+". Search for : "+Integer.toBinaryString(v)+" Pos "+ pos.length);
//         }else{
////            System.out.println(""+(i+1)+". Search for : "+Integer.toBinaryString(v)+" Not found");  
//         }
//         
//         }
//         System.out.println("Found for : "+found); 
         
//         int a= -100;
//         long b= ((long)a<<32)+3;//(a<<32);
//        
//         System.out.println("A = "+Integer.toBinaryString(a));
//         System.out.println("B = "+Long.toBinaryString(b));
         
         
         
         
    }
}
