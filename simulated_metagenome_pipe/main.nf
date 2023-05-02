process fetch_refs {
    conda "$projectDir/environment.yml"

    publishDir "${params.out_dir}/references/", mode: 'copy'

    maxForks 2

    input:
        val(row)

    output:
        tuple val(row), path("*.fasta"), stdout

    script:
        """
        ref_fetch.py --accession ${row.ref_accession} --email ${params.email} --total_reads ${params.total_reads} --proportion ${row.proportion}
        """
}

process gen_reads {
    conda "$projectDir/environment.yml"
    input:
        tuple val(row), path(ref_fasta), val(n_reads)
    
    output:
        tuple val(row), path(ref_fasta), emit: tax_metadata
        path("reads.fastq"), emit: tax_fastq

    script:

        """
        badread simulate --reference "${ref_fasta}" --quantity ${n_reads} > "reads.fastq"
        """

}


workflow {
    Channel
        .fromPath("${params.meta_manifest}")
        .splitCsv(header: true, sep: "\t")
        .set {manifest_ch}
    
    fetch_refs(manifest_ch) | gen_reads

    gen_reads.out.tax_fastq
        .collectFile(name: "simulated_metagenome.fastq", newLine: true, storeDir: "${params.out_dir}")
}