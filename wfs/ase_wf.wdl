version 1.0

import "gatk_select_variants.wdl" as gatk_select_variants_task
import "star_align_reads.wdl" as star_align_reads_task
import "samtools_vw_filter.wdl" as samtools_vw_filter_task
import "gatk_asereadcount.wdl" as gatk_asereadcount_task

workflow ASE_pipeline {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "# Pipeline for Allele-specific expression (ASE) table creation"
    }
    parameter_meta {
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
        complete_vcf: {
            description: "Zipped VCF with variants to filter only SNVs",
            extension: ".vcf.gz"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        fastq_r1: {
            description: "Sample R1 FASTQ file",
            extension: ".R1.fastq.gz"
        }
        fastq_r2: {
            description: "Sample R2 FASTQ file",
            extension: ".R2.fastq.gz"
        }
        star_index: {
            description: "TAR zip reference files from Gencode",
            extension: ".tar.gz"
        }
        gatk_docker: {
            description: "GATK docker image (file-GVgyqvQ02k8XGzgY2V8BVVBG)"
        }
        star_docker: {
            description: "TAR zip docker image from STAR stored on DNAnexus (file-GVjJgfQ02k8bYxZ9z70g869k)",
            extension: ".tar.gz"
        }
        samtools_docker: {
            description: "TAR zip docker image from SAMtools stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        # gatk_select_variants
        File fasta
        File fasta_fai
        File fasta_dict
        File complete_vcf
        String prefix
        String? gatk_docker

        # star_align_reads
        File fastq_r1
        File fastq_r2
        File star_index
        String? star_docker
        # STAR options
        Int? runThreadN
        Int? sjdbOverhang
        Int? limitBAMsortRAM
        Int? outBAMsortingThreadN

        # samtools_vw_filter
        String? samtools_docker
    }

    call gatk_select_variants_task.select_variants as select_variants {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            complete_vcf=complete_vcf,
            prefix=prefix,
            gatk_docker=gatk_docker
    }
    call star_align_reads_task.aligh_reads as aligh_reads {
        input:
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            prefix=prefix,
            star_index=star_index,
            snv_vcf=select_variants.snv_vcf,
            star_docker=star_docker
    }
    call samtools_vw_filter_task.samtools_vw_filter as samtools_vw_filter {
        input:
            star_bam=aligh_reads.star_bam,
            prefix=prefix,
            samtools_docker=samtools_docker
    }
    call gatk_asereadcount_task.gatk_asereadcount as gatk_asereadcount {
        input:
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            vw_bam=samtools_vw_filter.vw_bam,
            vw_bai=samtools_vw_filter.vw_bai,
            snv_vcf=select_variants.snv_vcf,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    output {
        # select_variants
        File snv_vcf=select_variants.snv_vcf
        # aligh_reads
        File star_bam = aligh_reads.star_bam
        File sjdb_txt = aligh_reads.sjdb_txt
        File sjdb_tab = aligh_reads.sjdb_tab
        File read_counts = aligh_reads.read_counts
        File junctions = aligh_reads.junctions
        File junctions_pass1 = aligh_reads.junctions_pass1
        File junctions_pass1_log = aligh_reads.junctions_pass1_log
        Array[File] logs = aligh_reads.logs
        # samtools_vw_filter
        File vw_bam = samtools_vw_filter.vw_bam
        File vw_bai = samtools_vw_filter.vw_bai
        # gatk_asereadcount
        File ase_table = gatk_asereadcount.ase_table
    }
}
