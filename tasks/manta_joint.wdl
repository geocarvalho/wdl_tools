version 1.0

task manta_joint {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads."
    }
    input {
        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai
        File fasta
        File fasta_fai
        String prefix
        File? manta_docker
    }
    String actual_manta_docker=select_first([manta_docker, "quay.io/biocontainers/manta:1.6.0--h9ee0642_2"])
    command {
        set -uexo pipefail
        mkdir output/
        if [[ ~{fasta} == *gz ]]
        then
            reference=$(echo ~{fasta} | sed 's/.gz//g');
            gunzip ~{fasta};
        else
            reference=~{fasta}; fi
        configManta.py \
            --bam ~{proband_bam} \
            --bam ~{mother_bam} \
            --bam ~{father_bam} \
            --referenceFasta $reference \
            --runDir output/
        python output/runWorkflow.py
        mv output/results/stats/alignmentStatsSummary.txt ~{prefix}_alignmentStatsSummary.txt
        mv output/results/stats/svCandidateGenerationStats.tsv ~{prefix}_svCandidateGenerationStats.tsv
        mv output/results/stats/svCandidateGenerationStats.xml ~{prefix}_svCandidateGenerationStats.xml
        mv output/results/stats/svLocusGraphStats.tsv ~{prefix}_svLocusGraphStats.tsv
        mv output/results/variants/candidateSV.vcf.gz ~{prefix}_candidateSV.vcf.gz
        mv output/results/variants/candidateSV.vcf.gz.tbi ~{prefix}_candidateSV.vcf.gz.tbi
        mv output/results/variants/candidateSmallIndels.vcf.gz ~{prefix}_candidateSmallIndels.vcf.gz
        mv output/results/variants/candidateSmallIndels.vcf.gz.tbi ~{prefix}_candidateSmallIndels.vcf.gz.tbi
        mv output/results/variants/diploidSV.vcf.gz ~{prefix}_diploidSV.vcf.gz
        mv output/results/variants/diploidSV.vcf.gz.tbi ~{prefix}_diploidSV.vcf.gz.tbi
    }
    runtime {
        docker: actual_manta_docker
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
        File alignmentStatsSummary = "${prefix}_alignmentStatsSummary.txt"
        File svCandidateGenerationStats = "${prefix}_svCandidateGenerationStats.tsv"
        File svCandidateGenerationStats_xml = "${prefix}_svCandidateGenerationStats.xml"
        File svLocusGraphStats = "${prefix}_svLocusGraphStats.tsv"
        File candidateSV = "${prefix}_candidateSV.vcf.gz"
        File candidateSV_tbi = "${prefix}_candidateSV.vcf.gz.tbi"
        File candidateSmallIndels = "${prefix}_candidateSmallIndels.vcf.gz"
        File candidateSmallIndels_tbi = "${prefix}_candidateSmallIndels.vcf.gz.tbi"
        File diploidSV = "${prefix}_diploidSV.vcf.gz"
        File diploidSV_tbi = "${prefix}_diploidSV.vcf.gz.tbi"
    }
}

workflow manta_joint_wf {
    parameter_meta {
        fasta: {
            description: "Unziped reference FASTA file stored by DNAnexus",
            extension: ".fasta"
        }
        fasta_fai: {
            description: "Index for reference FASTA file stored by DNAnexus",
            extension: ".fai"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        proband_bam: {
            description: "Last proband's BAM file created by the pipeline",
            extension: ".bam"
        }
        proband_bai: {
            description: "Index for the proband's BAM file used",
            extension: ".bai"
        }
        mother_bam: {
            description: "Last mother's BAM file created by the pipeline",
            extension: ".bam"
        }
        mother_bai: {
            description: "Index for the mother's BAM file used",
            extension: ".bai"
        }
        father_bam: {
            description: "Last father's BAM file created by the pipeline",
            extension: ".bam"
        }
        father_bai: {
            description: "Index for the father's BAM file used",
            extension: ".bai"
        }
        manta_docker: {
            description: "TAR zip docker image from Picard stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai
        File fasta
        File fasta_fai
        String prefix
        File? manta_docker
    }
    call manta_joint {
        input:
            proband_bam = proband_bam,
            proband_bai = proband_bai,
            mother_bam = mother_bam,
            mother_bai = mother_bai,
            father_bam = father_bam,
            father_bai = father_bai,
            prefix = prefix,
            fasta = fasta,
            fasta_fai = fasta_fai,
            manta_docker = manta_docker
    }
    output {
        File alignmentStatsSummary = manta_joint.alignmentStatsSummary
        File svCandidateGenerationStats = manta_joint.svCandidateGenerationStats
        File svCandidateGenerationStats_xml = manta_joint.svCandidateGenerationStats_xml
        File svLocusGraphStats = manta_joint.svLocusGraphStats
        File candidateSV = manta_joint.candidateSV
        File candidateSV_tbi = manta_joint.candidateSV_tbi
        File candidateSmallIndels = manta_joint.candidateSmallIndels
        File candidateSmallIndels_tbi = manta_joint.candidateSmallIndels_tbi
        File diploidSV = manta_joint.diploidSV
        File diploidSV_tbi = manta_joint.diploidSV_tbi
    }
}
