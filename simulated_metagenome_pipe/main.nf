process fetch_refs {
    input:
        tuple val(taxon_name), val(ref_accession), val(proportion)

    output:
        tuple val(taxon_name), val(ref_accession), val(proportion), path("${taxon_name}.fasta")

    script:
        """
        ref_fetch.py --accession ${ref_accession} --email ${params.email} > ${taxon_name}.fasta
        """
}

process gen_reads {
    input:
        tuple val(taxon_name), val(ref_accession), val(proportion), path(ref_fasta)
    
    output:
        tuple val(taxon_name), val(ref_accession), val(proportion), path(ref_fasta), emit: tax_metadata
        path("${taxon_name}_reads.fastq"), emit: tax_fastq

    script:
        n_reads = proportion * params.total_reads

        """
        badread simulate --reference ${ref_fasta} --quantity ${n_reads} > ${taxon_name}_reads.fastq
        """

}


workflow {
    Channel
        .fromPath("${params.meta_manifest}")
        .splitCsv(header: true, sep: "\t")
        .set {manifest_ch}
    
    fetch_refs(row.taxon_name, row.ref_accession, row.proportion) | gen_reads

    gen_reads.out.tax_fastq
        .collectFile(name: "simulated_metagenome.fastq", newLine: true, storeDir: "${params.out_dir}")
}