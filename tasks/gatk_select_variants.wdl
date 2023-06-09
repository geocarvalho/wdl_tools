version 1.0

task select_variants {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK SelectVariants select only SNVs"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File complete_vcf
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.4.0.0"])
    command {
        set -uexo pipefail
        tabix ~{complete_vcf}
        gatk SelectVariants \
            --java-options "-Xmx14G" \
            -R ~{fasta} \
            -V ~{complete_vcf} \
            --restrict-alleles-to BIALLELIC \
            --select-type-to-include SNP \
            -O ~{prefix}.selectvariants.vcf.gz
        gunzip ~{prefix}.selectvariants.vcf.gz
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
        dx_timeout: "1H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File snv_vcf = "${prefix}.selectvariants.vcf"
    }
}

workflow select_variants_wf {
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
        gatk_docker: {
            description: "GATK docker image (file-GVgyqvQ02k8XGzgY2V8BVVBG)"
        }
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File complete_vcf
        String prefix
        String? gatk_docker
    }
    call select_variants {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            complete_vcf=complete_vcf,
            prefix=prefix,
            gatk_docker=gatk_docker
    }
    output {
        File snv_vcf=select_variants.snv_vcf
    }
}
