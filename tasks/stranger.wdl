version 1.0

task stranger {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Annotates output files from ExpansionHunter with the pathologic implications of the repeat sizes."
    }
    input {
        File repeat_vcf
        String prefix
        File? stranger_docker
    }
    String actual_stranger_docker=select_first([stranger_docker, "gvcn/stranger:v0.8.1"])
    command {
        set -uexo pipefail
        stranger -f /bin/stranger/stranger/resources/variant_catalog_hg38.json \
        ~{repeat_vcf} > ~{prefix}_repeats.stranger.vcf
    }
    runtime {
        docker: actual_stranger_docker
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
        File stranger_vcf = "${prefix}_repeats.stranger.vcf"
    }
}

workflow stranger_wf {
    parameter_meta {
        repeat_vcf: {
            description: "ExpansionRepeat output VCF",
            extension: [".vcf", ".vcf.gz"]
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        stranger_docker: {
            description: "TAR zip docker image from Stranger stored on DNAnexus (file-GXj2gPQ02k8zP1f1k1XXF5Z2)",
            extension: ".tar.gz"
        }
    }
    input {
        File repeat_vcf
        String prefix
        File? stranger_docker
    }
    call stranger {
        input:
            repeat_vcf = repeat_vcf,
            prefix = prefix,
            stranger_docker = stranger_docker
    }
    output {
        File cnv_score_txt = stranger.stranger_vcf
    }
}
