# QC_piplinev1.0
The quality control pipline for sequencing reads(DNA and RNA library) of illumina,Pacbio (RS/RSII/Sequel/Sequel II),Nanopore platform

<b>Usage</b>: QC_pipline_v1.py [-h] [-c {illumina,pacbio,nanopore}] [-cr CORRECTION]  
		[-s SEQ_TYPE] [-omic OMIC_TYPE]  
		[-fastp_p FASTP_PARAMETERS] [-ccs_p CCS_PARAMETERS]  
		[-lima_p LIMA_PARAMETERS]  
		[-isoseq3_p ISOSEQ3_PARAMETERS] 
		[-lordec_p LORDEC_CORRECT_PARAMETERS]  
		[-nanoplot_p NANOPLOT_PARAMETERS]  
		inputs outputs  

<b>positional arguments: </b> 
  <b>inputs</b>                The input directiory(rawdata),the suffixes of short  
                        reads in the directory must be _1/2.fq or _1/2.fq.gz  
                        or _1/2.fastq or _1/2.fastq.gz (paired-end reads),or  
                        .fq/.fq.gz/.fastq/.fastq.gz (single-end reads);the  
                        suffixes of pacbio long reads must be  
                        .subreads.bam(Sequel platform) or  
                        .1.bax.h5,.2.bax.h5,.3.bax.h5(RS platform);the  
                        suffiexes of nanopore long reads must be .fastq.gz  
  <b>outputs</b>               The output directiory(cleandata and report)  

<b>optional arguments: </b> 
  <b>-h, --help</b>            show this help message and exit  
  <b>-c {illumina,pacbio,nanopore}, --choice {illumina,pacbio,nanopore}</b>  
                        Choose a sequencing machine,illumina for short  
                        reads,pacbio and nanopore for long reads,default is  
                        illumina  
  <b>-cr CORRECTION, --correction CORRECTION</b>  
                        The illumina reads for pacbio long reads  
                        correction(used in CLR mode,the self-correction is  
                        enough for CCS mode in normal condition),FASTA/Q  
                        file(s),e.g. -cr "reads1.fa,reads2.fq,reads3.fq.gz" or  
                        -cr "read.fa"  
  <b>-s SEQ_TYPE, --seq_type SEQ_TYPE</b>  
                        Paired end(PE) or singe end(SE) short reads in  
                        illumina sequencing platform,default is PE  
  <b>-omic OMIC_TYPE, --omic_type OMIC_TYPE</b>  
                        The omics type:DNA(genome) or  
                        RNA(transciptome),default is DNA  
  <b>-fastp_p FASTP_PARAMETERS, --fastp_parameters FASTP_PARAMETERS</b>  
                        The parameters for fastp softwares,e.g: -fastp_p "-w  
                        1", the parameters:--html and --json are defaultly set  
                        and named by the samples names to avoid overlapping in  
                        this pipline,please do not set again  
  <b>-ccs_p CCS_PARAMETERS, --ccs_parameters CCS_PARAMETERS</b>  
                        The parameters for SMRTlinks ccs tools,defalut  
                        parameters in DNA,"-noPolish --minPasses 1" in RNA  
  <b>-lima_p LIMA_PARAMETERS, --lima_parameters LIMA_PARAMETERS</b>  
                        The parameters for SMRTlinks lima tools,default  
                        parameters in DNA,"--isoseq --no-pbi" in RNA  
  <b>-isoseq3_p ISOSEQ3_PARAMETERS, --isoseq3_parameters ISOSEQ3_PARAMETERS</b>  
                        The parameters for SMRTlinks isoseq3  
                        tools(refine/cluster/polish),seperated by  
                        semicolon,e.g. -isoseq3_p "-j 1;--verbose;-r  
                        0.99",default is ";--verbose;"  
  <b>-lordec_p LORDEC_CORRECT_PARAMETERS, --lordec_correct_parameters LORDEC_CORRECT_PARAMETERS</b>  
                        The parameters for lordec_correct tools,"-k 19 -s 3"  
                        in DNA,default parameters in RNA  
  <b>-nanoplot_p NANOPLOT_PARAMETERS, --nanoplot_parameters NANOPLOT_PARAMETERS</b>  
                        The parameters for Nanoplot software  
