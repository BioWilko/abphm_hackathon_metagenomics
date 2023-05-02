process fetch_refs {
    conda "$projectDir/environment.yml"

    input:
        val(row)

    output:
        tuple val(row), path("${row.taxon_name}.fasta")

    script:
        """
        ref_fetch.py --accession ${row.ref_accession} --email ${params.email} > "${row.taxon_name}.fasta"
        """
}

process gen_reads {
    conda "$projectDir/environment.yml"
    input:
        tuple val(row), path(ref_fasta)
    
    output:
        tuple val(row), path(ref_fasta), emit: tax_metadata
        path("${row.taxon_name}_reads.fastq"), emit: tax_fastq

    script:

        n_reads = Float.parseFloat(row.proportion) * params.total_reads

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