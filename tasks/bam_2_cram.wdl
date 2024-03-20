version 1.0


task samtools_bam2cram {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "https://davetang.org/muse/2014/09/26/bam-to-cram/"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        String prefix
        String? fmt_option
        String? threads
        String? samtools_docker
    }
    String actual_samtools_docker=select_first([samtools_docker, "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"])
    String actual_threads=select_first([threads, "14"])
    String actual_fmt_option=select_first([fmt_option, "normal"])
    command {
        set -uexo pipefail
        ls ~{fasta_fai} ~{fasta_dict}
        if [[ ~{fasta} == *gz ]]
        then
            reference=$(echo ~{fasta} | sed 's/.gz//g');
            gunzip ~{fasta};
        else
            reference=~{fasta}; fi
        samtools view -@ ~{actual_threads} -T $reference -C --output-fmt-option ~{actual_fmt_option} -o ~{prefix}.~{actual_fmt_option}.cram ~{bam_file}
    }
    runtime {
        docker: actual_samtools_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
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
        File normal_cram = "${prefix}.${actual_fmt_option}.cram"
    }
}

task samtools_cram2bam {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "https://davetang.org/muse/2014/09/26/bam-to-cram/"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File cram
        String prefix
        String? fmt_option
        String? threads
        String? samtools_docker
    }
    String actual_samtools_docker=select_first([samtools_docker, "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"])
    String actual_threads=select_first([threads, "14"])
    String actual_fmt_option=select_first([fmt_option, "normal"])
    command {
        set -uexo pipefail
        ls ~{fasta_fai} ~{fasta_dict}
        if [[ ~{fasta} == *gz ]]
        then
            reference=$(echo ~{fasta} | sed 's/.gz//g');
            gunzip ~{fasta};
        else
            reference=~{fasta}; fi
        samtools view -@ ~{actual_threads} -T $reference -C --input-fmt-option ~{fmt_option} -o ~{prefix}.~{actual_fmt_option}.bam ~{cram}
    }
    runtime {
        docker: actual_samtools_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
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
        File unnormal_bam = "${prefix}_${fmt_option}.bam"
    }
}

task samtools_index {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for SAMtools index BAM files"
    }
    input {
        File aligned_file
        String? threads
        String? samtools_docker
    }
    String actual_samtools_docker=select_first([samtools_docker, "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"])
    String actual_threads=select_first([threads, "14"])
    command {
        set -uexo pipefail
        samtools index -@ ~{actual_threads} ~{aligned_file}
    }
    runtime {
        docker: actual_samtools_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "4H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File index = "${aligned_file}.crai"
    }
}

workflow samtools_cram_wf {
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
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        bam_file: {
            description: "Last BAM file created by the pipeline",
            extension: ".bam"
        }
        samtools_docker: {
            description: "TAR zip docker image from Picard stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        String prefix
        String? fmt_option
        String? threads
        String? samtools_docker
    }
    call samtools_bam2cram {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = bam_file,
            threads = threads,
            fmt_option = fmt_option,
            prefix = prefix,
            samtools_docker = samtools_docker
    }
    call samtools_index {
        input:
            aligned_file = samtools_bam2cram.normal_cram,
            threads = threads,
            samtools_docker = samtools_docker
    }
    output {
        File alignment = samtools_bam2cram.normal_cram
        File alignment_index = samtools_index.index
    }
}
