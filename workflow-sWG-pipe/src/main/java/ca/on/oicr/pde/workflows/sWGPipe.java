package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class sWGPipe extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;

    // Input Data
    private String tumorBam;
    private String outputFilenamePrefix;
    private String pon; // can be custom 

    //Tools
    private String samtools = "/oicr/local/analysis/sw//samtools/samtools-1.2//bin/samtools";
    private String java = "/oicr/local/analysis/sw/gatk/GenomeAnalysisTK-3.5-0";
    private String picardJar = "/.mounts/labs/PDE/Modules/sw/picard/2.12.1/picard.jar";
    private String gatkDir = "/oicr/local/analysis/sw/gatk/GenomeAnalysisTK-3.5-0";
    private String ichorScript = "/.mounts/labs/TGL/gsi/tools/ichorCNA/scripts/runIchorCNA.R";
    private String readCounter = "/.mounts/labs/TGL/gsi/tools/hmmcopy_utils/bin/readCounter";
    private String ichorExtDataPath = "/.mounts/labs/TGL/gsi/databases/ichorCNA_ext";
   

    //Memory allocation
    private Integer sWGMem = 16;
    private String javaMem = "-Xmx8g";

    //path to bin
    private String rScript = "/.mounts/labs/PDE/Modules/sw/R/R-3.5.1/bin/Rscript";
    private String rLib = "/.mounts/labs/gsiprojects/gsi/tools/R_libs";

    //ref Data
    private String refFasta = "/.mounts/labs/PDE/data/gatkAnnotationResources/hg19_random.fa";
    private String grch371000gIndels = "/.mounts/labs/TGL/gsi/databases/1000G_phase1.indels.b37.vcf_chr.vcf";
    private String grch37MillsIndels = "/.mounts/labs/TGL/gsi/databases/Mills_and_1000G_gold_standard.indels.b37.vcf_chr.vcf";
    private String dbsnpVcf = "/oicr/data/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP147/dbsnp147_chr.vcf";
    
    
    

    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(sWGPipe.class.getName());
    private String queue;

    // meta-types
    private final static String TXT_METATYPE = "text/plain";
    private final static String TAR_GZ_METATYPE = "application/tar-gzip";
    private final static String BAM_METATYPE = "application/bam";
    private final static String BAI_METATYPE = "application/bai";
    
    private String nt;
    private String windowSize = Integer.toString(1000000);
    private String maxCN = Integer.toString(5);
    private String includeHOMDFlag = "FALSE";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            tumorBam = getProperty("input_file_tumor");

            //Ext id
            outputFilenamePrefix = getProperty("external_name");
            
            //params
            nt = getProperty("num_threads");
            windowSize = getProperty("window_size");
            maxCN = getProperty("max_cn");
            includeHOMDFlag = getProperty("include_homd_flag");

            //samtools
            samtools = getProperty("samtools");
            java = getProperty("java");
            picardJar=getProperty("picard_jar");
            gatkDir=getProperty("gatk_dir");
            ichorScript=getProperty("ichorCNA_script");
            readCounter=getProperty("hmmcopy_readcounter");
            ichorExtDataPath=getProperty("ichorCNA_external_data");

            // ref
            refFasta = getProperty("ref_fasta");
            grch371000gIndels = getProperty("kg_indels");
            grch37MillsIndels = getProperty("mills_indels");
            dbsnpVcf=getProperty("dnsnp_vcf");
            pon = getProperty("panel_of_normals");

            //r path
            rScript = getProperty("rpath") + "/bin/Rscript";
            rLib = getProperty("rLib");
            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            sWGMem = Integer.parseInt(getProperty("swg_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("tumor");
        file0.setSourcePath(tumorBam);
        file0.setType("application/bam");
        file0.setIsInput(true);
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {


        // workflow : read inputs tumor sWG bam; run bam preprcoessing; pfo processed bam and dupmetrics
        // run ichorCNA;  pfo seg file, params file and zipped folder of all results
        Job parentJob = null;
        
        // copy to tmp dir 
        Job accessData = getWorkflow().createBashJob("access_bam");
        Command cmd = accessData.getCommand();
        cmd.addArgument("cp " + getFiles().get("tumor").getProvisionedPath() + " " + this.dataDir + this.outputFilenamePrefix + ".bam");
        accessData.setMaxMemory(Integer.toString(sWGMem * 1024));
        accessData.setQueue(getOptionalProperty("queue", ""));
        parentJob = accessData;
        // index bam
        String inputBam = this.dataDir + this.outputFilenamePrefix + ".bam";
        Job indexInputBam = this.bamIndex(inputBam);
        indexInputBam.addParent(parentJob);
        parentJob = indexInputBam;
        // pre processing 
        Job markDuplicates = this.markDuplicates(inputBam);
        markDuplicates.addParent(parentJob);
        parentJob = markDuplicates;
        
        // pfo 
        String dupMetrics = this.dataDir + this.outputFilenamePrefix + ".dup_metrics.txt";
        SqwFile dupMetricsTxt = createOutputFile(dupMetrics, TXT_METATYPE, this.manualOutput);
        dupMetricsTxt.getAnnotations().put("processed_bam_file", "sWG_BamPreProcess");
        markDuplicates.addFile(dupMetricsTxt);
    
        String dedupBam = this.dataDir + this.outputFilenamePrefix + ".dedup.bam";
        
        Job indexDedupBam = this.bamIndex(dedupBam);
        indexDedupBam.addParent(parentJob);
        parentJob = indexDedupBam;
        
        Job indelRealignerTargetCreator = this.indelRealignerTargetCreator(dedupBam);
        indelRealignerTargetCreator.addParent(parentJob);
        parentJob = indelRealignerTargetCreator;
        
        Job indelRealigner = this.indelRealigner(dedupBam);
        indelRealigner.addParent(parentJob);
        parentJob = indelRealigner;
        
        String dedupIndelRealignBam = this.dataDir + this.outputFilenamePrefix + ".dedup.realign.bam";
        Job baseRecalibrator = this.baseRecalibrator(dedupIndelRealignBam);
        baseRecalibrator.addParent(parentJob);
        parentJob = baseRecalibrator;
        
        Job bqsr = this.printReads(dedupIndelRealignBam);
        bqsr.addParent(parentJob);
        parentJob = bqsr;
        
        String dedupRealignRecalBam = this.dataDir + this.outputFilenamePrefix + ".dedup.realign.recal.bam"; // PFO
        SqwFile dedupRealignRecalBamFile = createOutputFile(dedupRealignRecalBam, BAM_METATYPE, this.manualOutput);
        dedupRealignRecalBamFile.getAnnotations().put("processed_bam_file", "sWG_BamPreProcess");
        bqsr.addFile(dedupRealignRecalBamFile);
        
        String dedupRealignRecalBai = this.dataDir + this.outputFilenamePrefix + ".dedup.realign.recal.bai"; // PFO
        SqwFile dedupRealignRecalBaiFile = createOutputFile(dedupRealignRecalBai, BAI_METATYPE, this.manualOutput);
        dedupRealignRecalBaiFile.getAnnotations().put("index_processed_bam_file", "sWG_BamPreProcess");
        bqsr.addFile(dedupRealignRecalBaiFile);
        
        // collect picard metrics 
        Job collectWGSmetrics = this.collectPicardMetrics(dedupRealignRecalBam);
        collectWGSmetrics.addParent(parentJob);
        
        Job collectInsertSizeMetrics = this.collectInsertSizeMetrics(dedupRealignRecalBam);
        collectInsertSizeMetrics.addParent(parentJob);
        
        // ichorCNA
        
        Job generateWigFile = this.wigGenerator(dedupRealignRecalBam);
        generateWigFile.addParent(parentJob);
        
        String wigFile = this.dataDir + this.outputFilenamePrefix + ".wig";
        Job runIchorCNA = this.runIchorCNA(wigFile);
        runIchorCNA.addParent(generateWigFile);
        
        Job zipOutp = this.zipIchorCNAOutp(this.dataDir + this.outputFilenamePrefix);
        zipOutp.addParent(runIchorCNA);
        
        // PFO
        String segFile = this.dataDir + this.outputFilenamePrefix + ".seg";
        SqwFile segFileSqw = createOutputFile(segFile, TXT_METATYPE, this.manualOutput);
        segFileSqw.getAnnotations().put("segment_data_from_the_tool", "ichorCNA");
        runIchorCNA.addFile(segFileSqw);
        
        String cnSegFile = this.dataDir + this.outputFilenamePrefix + ".cna.seg";
        SqwFile cnSegFileSqw = createOutputFile(cnSegFile, TXT_METATYPE, this.manualOutput);
        cnSegFileSqw.getAnnotations().put("cna_segment_data_from_the_tool", "ichorCNA");
        runIchorCNA.addFile(cnSegFileSqw);
        
        String segText = this.dataDir + this.outputFilenamePrefix + ".seg.txt";
        SqwFile segTextSqw = createOutputFile(segText, TXT_METATYPE, this.manualOutput);
        segTextSqw.getAnnotations().put("segment_data_text_from_the_tool", "ichorCNA");
        runIchorCNA.addFile(segTextSqw);
        
        String paramsText = this.dataDir + this.outputFilenamePrefix + ".params.txt";
        SqwFile paramsTextSqw = createOutputFile(paramsText, TXT_METATYPE, this.manualOutput);
        paramsTextSqw.getAnnotations().put("params_text_from_the_tool", "ichorCNA");
        runIchorCNA.addFile(paramsTextSqw);
        
        String correctedDepth = this.dataDir + this.outputFilenamePrefix + ".correctedDepth.txt";
        SqwFile correctedDepthSqw = createOutputFile(correctedDepth, TXT_METATYPE, this.manualOutput);
        correctedDepthSqw.getAnnotations().put("corrected_depth_from_the_tool", "ichorCNA");
        runIchorCNA.addFile(correctedDepthSqw);
        
        String tarFileOutp = this.dataDir + this.outputFilenamePrefix + ".tar.gz";
        SqwFile tarFileOutpSqw = createOutputFile(tarFileOutp, TAR_GZ_METATYPE, this.manualOutput);
        tarFileOutpSqw.getAnnotations().put("CNV_plots_folder", "ichorCNA");
        zipOutp.addFile(tarFileOutpSqw);
        
    }

    private Job markDuplicates(String inputBam) {
        Job markDuplicates = getWorkflow().createBashJob("mark_duplicates");
        Command cmd = markDuplicates.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar");
        cmd.addArgument(this.picardJar);
        cmd.addArgument("MarkDuplicates AS=TRUE");
        cmd.addArgument("I=" + inputBam);
        cmd.addArgument("M=" + this.dataDir + this.outputFilenamePrefix+".dup_metrics.txt");
        cmd.addArgument("O=" + this.dataDir + this.outputFilenamePrefix+".dup.bam");
        markDuplicates.setMaxMemory(Integer.toString(sWGMem * 1024));
        markDuplicates.setQueue(getOptionalProperty("queue", ""));
        return markDuplicates;
    }

    private Job bamIndex(String bamFile) {
        Job indexBam = getWorkflow().createBashJob("index_bam");
        Command cmd = indexBam.getCommand();
        cmd.addArgument(this.samtools + " index " + bamFile);
        indexBam.setMaxMemory(Integer.toString(sWGMem * 1024));
        indexBam.setQueue(getOptionalProperty("queue", ""));
        return indexBam;
    }

    private Job indelRealignerTargetCreator(String dedupBam) {
        Job indelRealignerTargetCreator = getWorkflow().createBashJob("realigner_target_creator");
        Command cmd = indelRealignerTargetCreator.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-Djava.io.tmpdir=" + this.tmpDir);
        cmd.addArgument("-jar");
        cmd.addArgument(this.gatkDir + "/GenomeAnalysisTK.jar");
        cmd.addArgument("-T RealignerTargetCreator -nt " + this.nt ); // set default nt = 8
        cmd.addArgument("-I " + dedupBam);
        cmd.addArgument("-R " +  this.refFasta);
        cmd.addArgument("-o " + this.dataDir + this.outputFilenamePrefix + ".intervals");
        cmd.addArgument("-dt NONE");
        cmd.addArgument("-known " + this.grch371000gIndels);
        cmd.addArgument("-known " + this.grch37MillsIndels);
        indelRealignerTargetCreator.setMaxMemory(Integer.toString(sWGMem * 1024));
        indelRealignerTargetCreator.setQueue(getOptionalProperty("queue", ""));
        return indelRealignerTargetCreator;
    }

    private Job indelRealigner(String dedupBam) {
        Job indelRealigner = getWorkflow().createBashJob("indel_realignment");
        Command cmd = indelRealigner.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-Djava.io.tmpdir=" + this.tmpDir);
        cmd.addArgument("-jar");
        cmd.addArgument(this.gatkDir + "/GenomeAnalysisTK.jar");
        cmd.addArgument("-T IndelRealigner");
        cmd.addArgument("-R " + this.refFasta);
        cmd.addArgument("-I " + dedupBam);
        cmd.addArgument("-targetIntervals " + this.dataDir +this.outputFilenamePrefix +".intervals");
        cmd.addArgument("-dt NONE");
        cmd.addArgument("-o " + this.dataDir + this.outputFilenamePrefix + ".dedup.realign.bam");
        cmd.addArgument("-known " + this.grch371000gIndels);
        cmd.addArgument("-known " + this.grch37MillsIndels);
        indelRealigner.setMaxMemory(Integer.toString(sWGMem * 1024));
        indelRealigner.setQueue(getOptionalProperty("queue", ""));
        return indelRealigner;
    }

    private Job baseRecalibrator(String dedupeRealignBam) {
        Job baseRecalibrator = getWorkflow().createBashJob("base_recalibration");
        Command cmd = baseRecalibrator.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-Djava.io.tmpdir=" + this.tmpDir);
        cmd.addArgument("-jar");
        cmd.addArgument(this.gatkDir + "/GenomeAnalysisTK.jar");
        cmd.addArgument("-T BaseRecalibrator -nct" + this.nt);
        cmd.addArgument("-I " + dedupeRealignBam );
        cmd.addArgument("-R " + this.refFasta);
        cmd.addArgument("-o " + this.dataDir + this.outputFilenamePrefix + ".dedup.realign.recal.grp");
        cmd.addArgument("-knownSites " + this.dbsnpVcf);
        cmd.addArgument("-rf BadCigar");
        cmd.addArgument("-cov ReadGroupCovariate");
        cmd.addArgument("-cov ContextCovariate");
        cmd.addArgument("-cov CycleCovariate");
        cmd.addArgument("-cov QualityScoreCovariate");
        cmd.addArgument("-dt NONE");
        cmd.addArgument("-knownSites " + this.grch371000gIndels);
        cmd.addArgument("-knownSites " + this.grch37MillsIndels);
        baseRecalibrator.setMaxMemory(Integer.toString(sWGMem * 1024));
        baseRecalibrator.setQueue(getOptionalProperty("queue", ""));
        return baseRecalibrator;
    }

    private Job printReads(String dedupeRealignBam) {
        Job printReads = getWorkflow().createBashJob("recaliberate_bam");
        Command cmd = printReads.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-Djava.io.tmpdir=" + this.tmpDir);
        cmd.addArgument("-jar");
        cmd.addArgument(this.gatkDir + "/GenomeAnalysisTK.jar");
        cmd.addArgument("-T PrintReads -nct " + this.nt);
        cmd.addArgument("-I " + dedupeRealignBam );
        cmd.addArgument("-R " + this.refFasta);
        cmd.addArgument("-BQSR " + this.dataDir + this.outputFilenamePrefix +  ".dedup.realign.recal.bam");
        cmd.addArgument("-rf BadCigar");
        cmd.addArgument("-dt NONE");
        printReads.setMaxMemory(Integer.toString(sWGMem * 1024));
        printReads.setQueue(getOptionalProperty("queue", ""));
        return printReads;
    }
    
    // picard tools
    private Job collectPicardMetrics(String dedupeRealignRecalBam) {
        Job collectPicardMetrics = getWorkflow().createBashJob("collect_wgs_metrics");
        Command cmd = collectPicardMetrics.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar");
        cmd.addArgument(this.picardJar);
        cmd.addArgument("CollectWgsMetrics");
        cmd.addArgument("I=" + dedupeRealignRecalBam );
        cmd.addArgument("O=" + this.dataDir + this.outputFilenamePrefix + ".wgs_metrics.txt");
        collectPicardMetrics.setMaxMemory(Integer.toString(sWGMem * 1024));
        collectPicardMetrics.setQueue(getOptionalProperty("queue", ""));
        return collectPicardMetrics;
    }
    
    private Job collectInsertSizeMetrics(String dedupeRealignRecalBam) {
        Job collectInsertSizeMetrics = getWorkflow().createBashJob("collect_insertsize_metrics");
        Command cmd = collectInsertSizeMetrics.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar");
        cmd.addArgument(this.picardJar);
        cmd.addArgument("CollectInsertSizeMetrics");
        cmd.addArgument("I=" + dedupeRealignRecalBam );
        cmd.addArgument("O=" + this.dataDir + this.outputFilenamePrefix + ".isize_metrics.txt");
        collectInsertSizeMetrics.setMaxMemory(Integer.toString(sWGMem * 1024));
        collectInsertSizeMetrics.setQueue(getOptionalProperty("queue", ""));
        return collectInsertSizeMetrics;
    }
    
    // ichorCNA
    private Job wigGenerator(String dedupeRealignRecalBam) {
        Job wigGenerator = getWorkflow().createBashJob("generate_wig");
        Command cmd = wigGenerator.getCommand();
        cmd.addArgument(this.readCounter);
        cmd.addArgument("--window " + this.windowSize); //default = 1000000 (bases)
        cmd.addArgument("-chromosome \"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY\"" );
        cmd.addArgument(dedupeRealignRecalBam);
        cmd.addArgument("|sed s/chrom=chr/chrom=/ >" + this.dataDir + this.outputFilenamePrefix + ".wig");
        wigGenerator.setMaxMemory(Integer.toString(sWGMem * 1024));
        wigGenerator.setQueue(getOptionalProperty("queue", ""));
        return wigGenerator;
    }
    
    
    private Job runIchorCNA(String wigFile) {
        Job runIchorCNA = getWorkflow().createBashJob("run_ichorCNA");
        Command cmd = runIchorCNA.getCommand();
        cmd.addArgument("R_LIBS=" + this.rLib);
        cmd.addArgument(this.rScript);
        cmd.addArgument(this.ichorScript);
        cmd.addArgument("--id " + this.outputFilenamePrefix);
        cmd.addArgument("--WIG " + wigFile);
        cmd.addArgument("--ploidy \"c(2,3)\"");
        cmd.addArgument("--normal \"c(0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9)\"");
        cmd.addArgument("--maxCN " + this.maxCN); // default = 5
        cmd.addArgument("--gcWig " + this.ichorExtDataPath + "/gc_hg19_1000kb.wig");
        cmd.addArgument("--mapWig " + this.ichorExtDataPath + "/map_hg19_1000kb.wig");
        cmd.addArgument("--centromere " + this.ichorExtDataPath + "/GRCh37.p13_centromere_UCSC-gapTable.txt");
        cmd.addArgument("--normalPanel " + this.pon);
        cmd.addArgument("--includeHOMD " + this.includeHOMDFlag); //defaul = FALSE
        cmd.addArgument("--chrs \"c(1:22, 'X')\"");
        cmd.addArgument("--chrTrain \"c(1:22)\"");
        cmd.addArgument("--estimateNormal True" );
        cmd.addArgument("--estimatePloidy True");
        cmd.addArgument("--estimateScPrevelance True");
        cmd.addArgument("--scStates \"c(1,3)\"");
        cmd.addArgument("--txnE 0.9999");
        cmd.addArgument("--txnStrength 10000");
        cmd.addArgument("--outDir " + this.dataDir + this.outputFilenamePrefix);
        runIchorCNA.setMaxMemory(Integer.toString(sWGMem * 1024));
        runIchorCNA.setQueue(getOptionalProperty("queue", ""));
        return runIchorCNA;
    }
    
    private Job zipIchorCNAOutp(String ichorResultsDir) {
        Job zipIchorCNAOutp = getWorkflow().createBashJob("zip_ichorCNA_plots");
        Command cmd = zipIchorCNAOutp.getCommand();
        cmd.addArgument("tar -czvf " +  ichorResultsDir + ".tar.gz" + " " + ichorResultsDir);
        zipIchorCNAOutp.setMaxMemory(Integer.toString(sWGMem * 1024));
        zipIchorCNAOutp.setQueue(getOptionalProperty("queue", ""));
        return zipIchorCNAOutp;
    }

}