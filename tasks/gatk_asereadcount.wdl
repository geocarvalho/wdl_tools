version 1.0

task gatk_asereadcount {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK ASEReadCounter table creation"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File vw_bam
        File vw_bai
        File snv_vcf
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.4.0.0"])
    command {
        set -uexo pipefail
        bgzip -c ~{snv_vcf} > ~{snv_vcf}.gz
        tabix -p vcf ~{snv_vcf}.gz
        gatk ASEReadCounter \
            --java-options "-Xmx14G" \
            -R ~{fasta} \
            -I ~{vw_bam} \
            -V ~{snv_vcf}.gz \
            -O ~{prefix}.wasp_filter.ASEReadCounter.table
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "5H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File ase_table = "${prefix}.wasp_filter.ASEReadCounter.table"
    }
}

workflow gatk_asereadcount_wf {
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
        vw_bam: {
            description: "BAM aligned by STAR alignReads with vW filter",
            extension: ".bam"
        }
        vw_bai: {
            description: "Index for BAM aligned by STAR alignReads with vW filter",
            extension: ".bai"
        }
        snv_vcf: {
            description: "VCF with only SNVs from same sample's WGS data",
            extension: ".vcf"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        gatk_docker: {
            description: "GATK docker image (file-GVgyqvQ02k8XGzgY2V8BVVBG)"
        }
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File vw_bam
        File vw_bai
        File snv_vcf
        String prefix
        String? gatk_docker
    }
    call gatk_asereadcount {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            vw_bam=vw_bam,
            vw_bai=vw_bai,
            snv_vcf=snv_vcf,
            prefix=prefix,
            gatk_docker=gatk_docker
    }
    output {
        File ase_table = gatk_asereadcount.ase_table
    }
}
