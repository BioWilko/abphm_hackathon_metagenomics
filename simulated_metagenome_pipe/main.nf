process fetch_refs {
    conda "$projectDir/environment.yml"

    input:
        val(row)

    output:
        tuple val(row), path("${row.ref_accession}_ref.fasta"), stdout

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
        path("${row.taxon_name}_reads.fastq"), emit: tax_fastq

    script:

        """
        badread simulate --reference "${ref_fasta}" --quantity ${n_reads} > "${row.taxon_name}_reads.fastq"
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