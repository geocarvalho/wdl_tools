version 1.0

task merge_fastqs {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task to run python script to merge FASTQ files"
    }
    parameter_meta {
        fastq_list: {
            description: "List of zipped FASTQ files",
            extension: ".fastq.gz"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        docker_image: {

        }
    }
    input {
        Array[File] fastq_list
        String prefix
        String? docker_image
    }
    String actual_docker_image=select_first([docker_image, "gvcn/merge_fastq:v0.0.4"])
    command {
        set -uexo pipefail
        python /home/bin/merge_fastqs.py -f ~{sep="," fastq_list} \
            -p ~{prefix}
    }
    runtime {
        docker: actual_docker_image
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
        File fastq_r1 = "${prefix}_R1_001_merge.fastq.gz"
        File fastq_r2 = "${prefix}_R2_001_merge.fastq.gz"
    }
}

workflow merge_fastq_wf {
    input {
        Array[File] fastq_list
        String prefix
        String? docker_image
    }
    call merge_fastqs {
        input:
            fastq_list=fastq_list,
            prefix=prefix,
            docker_image=docker_image
    }
    output {
        File fastq_r1=merge_fastqs.fastq_r1
        File fastq_r2=merge_fastqs.fastq_r2
    }
}
