version 1.0

import "../tasks/fastp.wdl" as fastp_task
import "../tasks/picard.wdl" as picard_task
import "../tasks/multiqc.wdl" as multiqc_task

workflow wgs_qc_wf {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Pipeline for WGS QC stats check."
    }
    parameter_meta {
        fastq_r1: {
            description: "Sample R1 FASTQ file",
            extension: ".R1.fastq.gz"
        }
        fastq_r2: {
            description: "Sample R2 FASTQ file",
            extension: ".R2.fastq.gz"
        }
        fastp_docker: {
            description: "TAR zip docker image from fastp stored on DNAnexus (file-GVzvBQ802k8k4zPYKqpfz1vZ)",
            extension: ".tar.gz"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        fasta: {
            description: "Unziped reference FASTA file stored by DNAnexus",
            extension: ".fasta"
        }
        fasta_fai: {
            description: "Index for reference FASTA file stored by DNAnexus",
            extension: ".fai"
        }
        fasta_dict: {
            description: "Dictionary for reference FASTA file stored by DNAnexus",
            extension: ".dict"
        }
        bam_file: {
            description: "Last BAM file created by the pipeline",
            extension: ".bam"
        }
        bam_index: {
            description: "Index for the BAM file used",
            extension: ".bai"
        }
        picard_docker: {
            description: "TAR zip docker image from Picard stored on DNAnexus ()",
            extension: ".tar.gz"
        }
        multiqc_image: {
            description: "TAR zip docker image from MultiQC stored on DNAnexus (file-GX16qZ802k8xKbkF797xKqP4)"
        }
        java_mem: {
            description: "String for Java memory (default: -Xmx100g)"
        }
    }
    input {
        # fastp
        File fastq_r1
        File fastq_r2
        String prefix
        String? fastp_docker
        # picard
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        File bam_index
        String? java_mem
        String? picard_docker
        String? multiqc_image
    }
    call fastp_task.fastp_stats as fastp_stats {
        input:
            fastq_r1 = fastq_r1,
            fastq_r2 = fastq_r2,
            prefix = prefix,
            fastp_docker = fastp_docker
    }
    call picard_task.picard_CollectMultipleMetrics as picard_CollectMultipleMetrics {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = bam_file,
            bam_index = bam_index,
            prefix = prefix,
            java_mem = java_mem,
            picard_docker = picard_docker
    }
    call picard_task.picard_CollectWgsMetrics as picard_CollectWgsMetrics {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = bam_file,
            bam_index = bam_index,
            prefix = prefix,
            java_mem = java_mem,
            picard_docker = picard_docker
    }
    call multiqc_task.multiqc as multiqc {
        input:
            stats_files = [
                fastp_stats.fastp_html, 
                fastp_stats.fastp_json, 
                picard_CollectMultipleMetrics.alignment_summary_metrics,
                picard_CollectMultipleMetrics.base_distribution_by_cycle,
                picard_CollectMultipleMetrics.base_distribution_by_cycle_metrics,
                picard_CollectMultipleMetrics.insert_size_histogram,
                picard_CollectMultipleMetrics.insert_size_metrics,
                picard_CollectMultipleMetrics.quality_by_cycle,
                picard_CollectMultipleMetrics.quality_by_cycle_metrics,
                picard_CollectMultipleMetrics.quality_distribution,
                picard_CollectMultipleMetrics.quality_distribution_metrics,
                picard_CollectMultipleMetrics.read_length_histogram,
                picard_CollectWgsMetrics.collect_wgs_metrics
            ],
            prefix = prefix,
            multiqc_image = multiqc_image
    }
    output {
        File fastp_html = fastp_stats.fastp_html
        File fastp_json = fastp_stats.fastp_json
        File alignment_summary_metrics = picard_CollectMultipleMetrics.alignment_summary_metrics
        File base_distribution_by_cycle = picard_CollectMultipleMetrics.base_distribution_by_cycle
        File base_distribution_by_cycle_metrics = picard_CollectMultipleMetrics.base_distribution_by_cycle_metrics
        File insert_size_histogram = picard_CollectMultipleMetrics.insert_size_histogram
        File insert_size_metrics = picard_CollectMultipleMetrics.insert_size_metrics
        File quality_by_cycle = picard_CollectMultipleMetrics.quality_by_cycle
        File quality_by_cycle_metrics = picard_CollectMultipleMetrics.quality_by_cycle_metrics
        File quality_distribution = picard_CollectMultipleMetrics.quality_distribution
        File quality_distribution_metrics = picard_CollectMultipleMetrics.quality_distribution_metrics
        File read_length_histogram = picard_CollectMultipleMetrics.read_length_histogram
        File collect_wgs_metrics = picard_CollectWgsMetrics.collect_wgs_metrics
        File multiqc_html = multiqc.multiqc_html
        File multiqc_data = multiqc.multiqc_data
    }
}
