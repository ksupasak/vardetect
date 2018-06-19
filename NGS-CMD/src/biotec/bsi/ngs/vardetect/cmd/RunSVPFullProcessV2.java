/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence;
import biotec.bsi.ngs.vardetect.core.VariationResult;
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
public class RunSVPFullProcessV2 {
    
    private static int minPeakF=5;
    private static int minPeakB=2;
    private static int merCov=50;
    private static int maxDup=1;
    private static int countDup=1;
    private static int numThread=1;
    private static long numSkip=0;
    private static int percentBaseRepeat=50;
    private static String refPath="";
    private static String gffFile="";
    private static String inputPath = "";
    private static String samtools = "";
    private static String mainBamFile = "";
    private static int numRead=1000000;
    private static int numMer=16;
    private static int filterMode=3;
    private static int minIndelSize=30;
    private static int minReadPerGroup=5;
    private static int maxExtend=100;
    
    public static void main(String[] args) throws IOException, FileNotFoundException, InterruptedException{
    // create the command line parser
        CommandLineParser parser = new DefaultParser();

        // create the Options
        Options options = new Options();
//        options.addOption( "a", "all", false, "do not hide entries starting with ." );
//        options.addOption( "A", "almost-all", false, "do not list implied . and .." );
//        options.addOption( "b", "escape", false, "print octal escapes for nongraphic "
        options.addOption( Option.builder("b").longOpt("bam-file")
                                        .desc("input file in bam format (.bam)")
                                        .required()
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("B").longOpt("main bam-file")
                                        .desc("Main bam file that include both map and unmap (.bam)")
                                        .required()
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("S").longOpt("samtools path")
                                        .desc("Full path of samtools ex. use/bin/samtools")
                                        .required()
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("r").longOpt("ref-file")
                                        .desc("Reference file path")
                                        .required()
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("a").longOpt("gff3-file")
                                        .desc("annotation gff3 file path")
                                        .hasArg()
                                        .argName("file")
                                        .build());
        options.addOption( Option.builder("t").longOpt("thread-num")
                                        .desc("number of Thread [1]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("n").longOpt("read-num")
                                        .desc("total number of read per thread [1000000]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("s").longOpt("skip-read")
                                        .desc("number of read to skip (use when debug) [0]" )
                                        .hasArg()
                                        .argName("long")
                                        .build());
        options.addOption( Option.builder("d").longOpt("max-dup")
                                        .desc("maximum duplicate pattern [1]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("p").longOpt("mix-pattern")
                                        .desc("minimum peak pattern [10,5]" )
                                        .hasArgs().numberOfArgs(2)
                                        .argName("int int")
                                        .build());
        options.addOption( Option.builder("m").longOpt("mer-size")
                                        .desc("size of mer [16]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("C").longOpt("mer-coverage")
                                        .desc("Filter value : number of mer coverage [50]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("F").longOpt("filter-mode")
                                        .desc("Filter mode has 1,2,3  number of mode [1]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("R").longOpt("ref-repeat")
                                        .desc("Filter value : percent of base repeat (percent of lowercase letter) [50]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("D").longOpt("filter-dup")
                                        .desc("Filter duplication count [1]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("i").longOpt("min-indelSize")
                                        .desc("minimum indel size [30]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("g").longOpt("min-readPerGroup")
                                        .desc("minimum read per group (clustering) [5]" )
                                        .hasArg()
                                        .argName("int")
                                        .build());
        options.addOption( Option.builder("e").longOpt("max-extendRef")
                                        .desc("maximum base extend per side [100]" )
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
            if(line.hasOption("a")) gffFile = line.getOptionValue("a");
            if(line.hasOption("b")) inputPath = line.getOptionValue("b");
            if(line.hasOption("R")) percentBaseRepeat = Integer.parseInt(line.getOptionValue("R"));
            if(line.hasOption("s")) numSkip = Long.parseLong(line.getOptionValue("s"));
            if(line.hasOption("t")) numThread = Integer.parseInt(line.getOptionValue("t"));
            if(line.hasOption("i")) minIndelSize = Integer.parseInt(line.getOptionValue("i"));
            if(line.hasOption("g")) minReadPerGroup = Integer.parseInt(line.getOptionValue("g"));
            if(line.hasOption("e")) maxExtend = Integer.parseInt(line.getOptionValue("e"));
            if(line.hasOption("B")) mainBamFile = line.getOptionValue("B");
            if(line.hasOption("S")) samtools = line.getOptionValue("S");

            if(line.hasOption("h")){
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp( "SVP parameter setting help", options );
            }
        }
        catch( ParseException exp ) {
            System.out.println( "Unexpected exception:" + exp.getMessage() );
        }
        
        
//        String refPath = args[0];
//         String refPath = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa";
//        String inBam = args[1];
        System.out.println("gffFile = "+gffFile);
        String[] dummy = inputPath.split("\\.bam");
        System.out.println("root="+dummy[0]);
        String outputAllFile = dummy[0]+"_all.out";
        String outputFilterFile = dummy[0]+"_filter.out";

        CombineReferenceSequence ref = SequenceUtil.getCombineReferenceSequence(refPath,numMer); //runFile hg19.fa

        ref.setMinimumPeakPattern(minPeakF, minPeakB);       // minimum mer per pattern A and B

// for large batch with multi thread         
        ref.setNumberOfThread(numThread);              // numthread
//         ref.setTotalRead(86000000);
//         ref.setTotalRead(8000);

        ref.setChunkRead(numRead);             // number of read per thread defualt 1 million
        ref.setSkipRead(numSkip);                    
        ref.setMaximumDuplicatePattern(maxDup);     // Consider repeat or not (1 is consider not repeat, 2 is consider repeat 2 location, 3 is consider repeat 3 location ... max is 9 location)

        ref.setFilterMerCoverage(merCov);          // defualt 50
        ref.setFilterMode(filterMode);                  // defualt 1
        ref.setFilterRefRepeatCount(percentBaseRepeat);       // defualt 50
        ref.setFilterDupCount(countDup);              // defualt 1;

        ref.setMinimumIndelSize(minIndelSize);            // defualt is 30 
// for large batch with multi thread         
//         ref.setNumberOfThread(4);
//         ref.setTotalRead(100000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(1);

// for small batch with multi thread         
//         ref.setNumberOfThread(4);
//         ref.setTotalRead(10000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(1);
//         ref.setRandomAccess(true);


// for debug         
//         ref.setNumberOfThread(1);
//         ref.setTotalRead(100);
//         ref.setSkipRead(10242);
//         ref.setSkipRead(51207);
////         ref.setSkipRead(86893);
//         ref.setSkipRead(32382);
//////         ref.setSkipRead(21555455);
//         ref.setMaximumDuplicatePattern(3);
//         ref.setMinimumPeakPattern(10, 2);
//         ref.setRandomAccess(true);


        ref.prepare();
        ref.setOutputFile(outputAllFile);
        ref.setOutputSVFile(outputFilterFile);
        ref.runProfileSV(inputPath);
        
        System.out.println("Post processing");
        String inputFile = outputFilterFile;
        String saveFile = outputFilterFile.split("\\.out")[0];
        String refFaIdx = refPath+".fai";
        System.out.println("refFaIdx = "+refFaIdx);
        VariationResult varRes = SequenceUtil.readVersion2AlignmentResult(inputFile);
        varRes.createRefIndex(refFaIdx);
        varRes.checkRefIndex();
        varRes.analyzeCoverage();
//        /**
//         * Rough SV protocol
//         */
//        varRes.classifyRoughSVType();
//        if(gffFile.equals("")){
//            varRes.writeStructureVariantV2SortedCoverageReportToFile(saveFile, minReadPerGroup);
//            varRes.writeStructureVariantV2SortedCoverageGroupInfoReportToFile(saveFile, minReadPerGroup);
//        }else{
//            varRes.writeStructureVariantV2SortedCoverageReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, minReadPerGroup, numMer);
//            varRes.writeStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, minReadPerGroup, numMer);
//        }
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refPath, refFaIdx, "TD", minReadPerGroup, maxExtend);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refPath, refFaIdx, "ID", minReadPerGroup, maxExtend);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refPath, refFaIdx, "IC", minReadPerGroup, maxExtend);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refPath, refFaIdx, "IT", minReadPerGroup, maxExtend);
//        /*****/
        
        /**
         * Precise SV protocol
         */
        varRes.classifyPreciseSVType(minReadPerGroup);
        varRes.identifyCorrectness(mainBamFile,samtools);
        varRes.sortInsertionByInsertJunctiontLowtoHigh();
        varRes.FilterIntraInsertion();  // Del Del filter
        
        if(!gffFile.equals("")){
            System.out.println("Generate Annotation report");
            varRes.writePreciseStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationExcel(inputFile, gffFile, refFaIdx, minReadPerGroup, numMer);
            varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refPath, refFaIdx, gffFile, "TD", minReadPerGroup, maxExtend, numMer);
            varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refPath, refFaIdx, gffFile, "D", minReadPerGroup, maxExtend, numMer);
            varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refPath, refFaIdx, gffFile, "IA", minReadPerGroup, maxExtend, numMer);
            varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refPath, refFaIdx, gffFile, "IE", minReadPerGroup, maxExtend, numMer);
            varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refPath, refFaIdx, gffFile, "CH", minReadPerGroup, maxExtend, numMer);
        }else{
            System.out.println("Generate report");
            varRes.writePreciseStructureVariantV2SortedCoverageGroupInfoReportExcel(inputFile, refFaIdx, minReadPerGroup);
            varRes.writeVisualizePreciseSVType(inputFile, refPath, refFaIdx, "TD", minReadPerGroup, maxExtend);
            varRes.writeVisualizePreciseSVType(inputFile, refPath, refFaIdx, "D", minReadPerGroup, maxExtend);
            varRes.writeVisualizePreciseSVType(inputFile, refPath, refFaIdx, "IA", minReadPerGroup, maxExtend);
            varRes.writeVisualizePreciseSVType(inputFile, refPath, refFaIdx, "IE", minReadPerGroup, maxExtend);
            varRes.writeVisualizePreciseSVType(inputFile, refPath, refFaIdx, "CH", minReadPerGroup, maxExtend);
        }
        
        /**********/
        
        System.out.println("Done");
        System.exit(0);
    }
}
