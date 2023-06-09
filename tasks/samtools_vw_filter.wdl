version 1.0

task samtools_vw_filter {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for SAMtools view filter vW tags"
    }
    input {
        File star_bam
        String prefix
        String? samtools_docker
    }
    String actual_samtools_docker=select_first([samtools_docker, "quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"])
    command {
        set -uexo pipefail
        samtools index -@ 30 ~{star_bam}
        samtools view -@ 30 ~{star_bam} -e '![vW] || [vW]==1' --with-header --bam --output ~{prefix}.Aligned.sortedByCoord.out.wasp_filter.bam
        samtools index -@ 30 ~{prefix}.Aligned.sortedByCoord.out.wasp_filter.bam
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
        dx_timeout: "2H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File vw_bam = "${prefix}.Aligned.sortedByCoord.out.wasp_filter.bam"
        File vw_bai = "${prefix}.Aligned.sortedByCoord.out.wasp_filter.bam.bai"
    }
}

workflow samtools_vw_filter_wf {
    parameter_meta {
        star_bam: {
            description: "BAM aligned by STAR alignReads",
            extension: ".bam"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        samtools_docker: {
            description: "TAR zip docker image from SAMtools stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        File star_bam
        String prefix
        String? samtools_docker
    }
    call samtools_vw_filter {
        input:
            star_bam=star_bam,
            prefix=prefix,
            samtools_docker=samtools_docker
    }
    output {
        File vw_bam = samtools_vw_filter.vw_bam
        File vw_bai = samtools_vw_filter.vw_bai
    }
}
