version 1.0

task run_bash_script {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task to run a script, normally used to bash scripts to download files."
    }
    parameter_meta {
        bash_script: {
            description: "Bash script you want to run.",
            extension: ".sh"
        }
    }
    input {
        File bash_script
    }
    command {
        set -uexo pipefail
        bash ~{bash_script}
    }
    runtime {
        dx_instance_type: "mem1_ssd1_v2_x8"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "9H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        Array[File] zip_files = glob("*.gz")
        Array[File] bam_files = glob("*.bam")
        Array[File] bai_files = glob("*.bai")
        Array[File] tbi_files = glob("*.tbi")
    }
}

workflow run_bash_script_wf {
    input {
        File bash_script
    }
    call run_bash_script {
        input:
            bash_script=bash_script
    }
    output {
        Array[File] zip_files = run_bash_script.zip_files
        Array[File] bam_files = run_bash_script.bam_files
        Array[File] bai_files = run_bash_script.bai_files
        Array[File] tbi_files = run_bash_script.tbi_files
    }
}
