version 1.0

task qualimap_bamqc {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Not tested "
    }
    input {
        File bam_file
        File bam_index
        String prefix
        String? memory
        String? threads
        String? qualimap_docker
    }
    String actual_memory=select_first([memory, "25G"])
    String actual_threads=select_first([threads, "14"])
    String actual_qualimap_docker=select_first([qualimap_docker, "gvcn/qualimap:2.2.2a"])
    command {
        set -uexo pipefail
        qualimap --java-mem-size=~{actual_memory} bamqc \
            -bam ~{bam_file} \
            -outdir ~{prefix}_qualimap \
            -nt ~{actual_threads}
        mv ~{prefix}_qualimap/genome_results.txt ~{prefix}_qualimap/~{prefix}_genome_results.txt
        mv ~{prefix}_qualimap/qualimapReport.html ~{prefix}_qualimap/~{prefix}_qualimapReport.html
    }
    runtime {
        docker: actual_qualimap_docker
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
        File genome_results = "${prefix}_qualimap/${prefix}_genome_results.txt"
        File qualimap_report = "${prefix}_qualimap/~{prefix}_qualimapReport.html"
    }
}

workflow qualimap_bamqc_wf {
    parameter_meta {
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        bam_file: {
            description: "Last BAM file created by the pipeline",
            extension: ".bam"
        }
        bam_index: {
            description: "Index for the BAM file used",
            extension: ".bai"
        }
        qualimap_docker: {
            description: "TAR zip docker image from Picard stored on DNAnexus ()",
            extension: ".tar.gz"
        }
        memory: {
            description: "Amount of memory to be used (default: 25G)"
        }
        threads: {
            description: "Number of threads to be used (default: 14)"
        }
    }
    input {
        File bam_file
        File bam_index
        String prefix
        String? memory
        String? threads
        String? qualimap_docker
    }
    call qualimap_bamqc {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            prefix = prefix,
            memory = memory,
            threads = threads,
            qualimap_docker = qualimap_docker
    }
    output {
        File genome_results = qualimap_bamqc.genome_results
        File qualimap_report = qualimap_bamqc.qualimap_report
    }
}
