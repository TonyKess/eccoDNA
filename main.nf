#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    eccoDNA - Nextflow DSL2
========================================================================================
    Author: Your Name
    Usage:
        nextflow run main.nf -c nextflow.config --input samplesheet.csv --marker FISHE
========================================================================================
*/

// Print pipeline header
def helpMessage() {
    log.info"""
    ========================================
    eccoDNA Metabarcoding Pipeline
    ========================================
    Usage:
        nextflow run main.nf -c nextflow.config \\
            --input samplesheet.csv \\
            --marker FISHE \\
            --outdir results
    
    Required arguments:
        --input         Path to sample sheet (CSV with columns: sample,read1,read2)
        --marker        Marker name (must match a marker defined in nextflow.config)
                        Currently configured: ${params.markers.keySet().join(', ')}
        --outdir        Output directory
    
    Optional arguments:
        --blast.database      Path to BLAST database (default: from config)
        --skip_blast          Skip BLAST step
        --skip_sss_filter     Skip automatic SSS filtering
        --sss.min_threshold   SSS threshold for filtering (default: 99.0)
        --help                Print this help message
    
    Marker-specific parameters (length ranges, primers) are defined in nextflow.config.
    To add a new marker, edit the 'markers' section in nextflow.config.
    ========================================
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Please provide --input samplesheet.csv"
}
if (!params.marker) {
    error "Please provide --marker"
}

// Get marker-specific parameters
def markerParams = params.markers[params.marker]
if (!markerParams) {
    error "Unknown marker: ${params.marker}. Available: ${params.markers.keySet()}"
}

log.info """
========================================
Pipeline Parameters
========================================
Project      : ${params.project_name}
Marker       : ${params.marker}
Input        : ${params.input}
Output dir   : ${params.outdir}
Length range : ${markerParams.min_length}-${markerParams.max_length} bp
BLAST DB     : ${params.blast.database}
========================================
"""

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process MERGE_READS {
    tag "$sample"
    publishDir "${params.outdir}/merged", mode: 'copy', pattern: "*.assembled.fastq"
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}.merged.assembled.fastq"), emit: merged
    path "${sample}.pear.log", emit: log
    
    script:
    """
    pear -f ${read1} -r ${read2} \\
         -o ${sample}.merged \\
         -q ${params.pear.quality_threshold} \\
         -v ${params.pear.min_overlap} \\
         -j ${task.cpus} \\
         > ${sample}.pear.log 2>&1
    """
}

process QUALITY_FILTER {
    tag "$sample"
    publishDir "${params.outdir}/filtered", mode: 'copy'
    
    input:
    tuple val(sample), path(merged)
    
    output:
    tuple val(sample), path("${sample}.filtered.fasta"), emit: filtered
    path "${sample}.filter.log", emit: log
    
    script:
    """
    vsearch --fastq_filter ${merged} \\
            --fastq_maxee ${params.vsearch.max_ee} \\
            --fastq_minlen ${params.vsearch.min_len} \\
            --fastq_maxlen ${params.vsearch.max_len} \\
            --relabel ${sample}. \\
            --sample ${sample} \\
            --fastaout ${sample}.filtered.fasta \\
            > ${sample}.filter.log 2>&1
    """
}

process CONCATENATE {
    publishDir "${params.outdir}/concatenated", mode: 'copy'
    
    input:
    path(filtered_files)
    
    output:
    path "all.filtered.fasta", emit: concatenated
    
    script:
    """
    cat ${filtered_files} > all.filtered.fasta
    """
}

process DEREPLICATE {
    publishDir "${params.outdir}/dereplicated", mode: 'copy'
    
    input:
    path(concatenated)
    
    output:
    path "all.derep.fasta", emit: dereplicated
    path "derep.log", emit: log
    
    script:
    """
    vsearch --derep_fulllength ${concatenated} \\
            --minuniquesize ${params.vsearch.minuniquesize} \\
            --sizeout \\
            --output all.derep.fasta \\
            > derep.log 2>&1
    """
}

process DENOISE {
    publishDir "${params.outdir}/denoised", mode: 'copy'
    
    input:
    path(dereplicated)
    
    output:
    path "unoise_centroids.fasta", emit: denoised
    path "denoise.log", emit: log
    
    script:
    """
    vsearch --cluster_unoise ${dereplicated} \\
            --minsize ${params.vsearch.minsize} \\
            --unoise_alpha ${params.vsearch.unoise_alpha} \\
            --centroids unoise_centroids.fasta \\
            --relabel ASV_ \\
            > denoise.log 2>&1
    """
}

process CHIMERA_REMOVAL {
    publishDir "${params.outdir}/chimera_free", mode: 'copy'
    
    input:
    path(denoised)
    
    output:
    path "ASVs_nochimera.fasta", emit: nochimera
    path "chimera.log", emit: log
    
    script:
    """
    vsearch --uchime3_denovo ${denoised} \\
            --nonchimeras ASVs_nochimera.fasta \\
            > chimera.log 2>&1
    """
}

process LENGTH_FILTER {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path(nochimera)
    val(marker)
    val(min_len)
    val(max_len)
    
    output:
    path "${marker}_ASVs.fasta", emit: asvs
    
    script:
    """
    seqkit seq ${nochimera} \\
           --min-len ${min_len} \\
           --max-len ${max_len} \\
           > ${marker}_ASVs.fasta
    
    # Report ASV count
    echo "Final ASV count: \$(grep -c '^>' ${marker}_ASVs.fasta)" > asv_count.txt
    """
}

process GENERATE_ASV_TABLE {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path(asvs)
    path(concatenated)
    val(marker)
    
    output:
    path "${marker}_ASV_table.tsv", emit: table
    path "${marker}_ASV_table.biom", emit: biom
    path "asv_table.log", emit: log
    
    script:
    """
    vsearch --usearch_global ${concatenated} \\
            --db ${asvs} \\
            --id ${params.vsearch.identity_threshold} \\
            --threads ${task.cpus} \\
            --otutabout ${marker}_ASV_table.tsv \\
            --biomout ${marker}_ASV_table.biom \\
            > asv_table.log 2>&1
    """
}

process BLAST_TAXONOMY {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path(asvs)
    val(marker)
    
    output:
    path "${marker}_blast_results.tsv", emit: blast
    path "blast.log", emit: log
    
    script:
    """
    export BLASTDB=${params.blast.database}
    blastn -db ${params.blast.db} \\
           -query ${asvs} \\
           -max_target_seqs ${params.blast.max_target_seqs} \\
           -out ${marker}_blast_results.tsv \\
           -evalue ${params.blast.evalue} \\
           -outfmt "${params.blast.outfmt}" \\
           -num_threads ${task.cpus} \\
           > blast.log 2>&1
    """
}

process FILTER_SSS {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path(blast_results)
    path(asv_table)
    val(marker)
    val(sss_threshold)
    
    output:
    path "${marker}_filtered_blast.tsv", emit: filtered_blast
    path "${marker}_asv_taxonomy.tsv", emit: taxonomy
    path "${marker}_reduced_asv_table.tsv", emit: reduced_table
    path "sss_filtering_report.txt", emit: report
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    
    # Read BLAST results
    blast_cols = ['qseqid', 'sseqid', 'pident', 'evalue', 'length', 'qcovs', 'ssciname', 'sblastname', 'scomname', 'staxid']
    blast_df = pd.read_csv('${blast_results}', sep='\\t', names=blast_cols, header=None)
    
    # Calculate SSS (Sequence Similarity Score)
    blast_df['sss'] = (blast_df['pident'] * blast_df['qcovs']) / 100
    
    # Filter by SSS threshold
    filtered_df = blast_df[blast_df['sss'] >= ${sss_threshold}].copy()
    
    # Save filtered BLAST results
    filtered_df.to_csv('${marker}_filtered_blast.tsv', sep='\\t', index=False)
    
    # Create taxonomy summary: ASV -> species hits with SSS >= threshold
    taxonomy_summary = filtered_df.groupby('qseqid').agg({
        'ssciname': lambda x: list(x.unique()),
        'sss': ['max', 'mean', 'count']
    }).reset_index()
    
    taxonomy_summary.columns = ['ASV', 'species_hits', 'max_sss', 'mean_sss', 'n_hits']
    taxonomy_summary['species_hits'] = taxonomy_summary['species_hits'].apply(lambda x: ';'.join(x))
    taxonomy_summary.to_csv('${marker}_asv_taxonomy.tsv', sep='\\t', index=False)
    
    # Read ASV table
    asv_table = pd.read_csv('${asv_table}', sep='\\t')
    asv_col = asv_table.columns[0]  # First column is ASV ID
    
    # Merge ASV table with taxonomy (only keep ASVs that pass SSS filter)
    reduced_table = asv_table[asv_table[asv_col].isin(taxonomy_summary['ASV'])]
    reduced_table = reduced_table.merge(taxonomy_summary[['ASV', 'n_hits', 'max_sss']], 
                                        left_on=asv_col, right_on='ASV', how='left')
    reduced_table.to_csv('${marker}_reduced_asv_table.tsv', sep='\\t', index=False)
    
    # Generate report
    n_total_asvs = len(asv_table)
    n_blast_hits = len(blast_df['qseqid'].unique())
    n_filtered_asvs = len(taxonomy_summary)
    pct_retained = (n_filtered_asvs / n_total_asvs * 100) if n_total_asvs > 0 else 0
    
    with open('sss_filtering_report.txt', 'w') as f:
        f.write(f"SSS Filtering Report - ${marker}\\n")
        f.write(f"{'='*50}\\n")
        f.write(f"SSS threshold: ${sss_threshold}%\\n")
        f.write(f"Total ASVs in table: {n_total_asvs}\\n")
        f.write(f"ASVs with BLAST hits: {n_blast_hits}\\n")
        f.write(f"ASVs passing SSS filter: {n_filtered_asvs}\\n")
        f.write(f"Retention rate: {pct_retained:.2f}%\\n")
        f.write(f"\\nTop species by number of ASVs:\\n")
        
        # Count ASVs per species
        species_counts = filtered_df.groupby('ssciname')['qseqid'].nunique().sort_values(ascending=False).head(20)
        for species, count in species_counts.items():
            f.write(f"  {species}: {count} ASVs\\n")
    
    print(f"SSS filtering complete:")
    print(f"  Total ASVs: {n_total_asvs}")
    print(f"  Retained after SSS >= ${sss_threshold}%: {n_filtered_asvs} ({pct_retained:.2f}%)")
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Read sample sheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            tuple(row.sample, file(row.read1), file(row.read2))
        }
        .set { samples_ch }
    
    // Step 1: Merge reads
    MERGE_READS(samples_ch)
    
    // Step 2: Quality filter
    QUALITY_FILTER(MERGE_READS.out.merged)
    
    // Step 3: Concatenate all filtered reads
    QUALITY_FILTER.out.filtered
        .map { sample, fasta -> fasta }
        .collect()
        .set { filtered_files_ch }
    
    CONCATENATE(filtered_files_ch)
    
    // Step 4: Dereplicate
    DEREPLICATE(CONCATENATE.out.concatenated)
    
    // Step 5: Denoise with UNOISE
    DENOISE(DEREPLICATE.out.dereplicated)
    
    // Step 6: Remove chimeras
    CHIMERA_REMOVAL(DENOISE.out.denoised)
    
    // Step 7: Length filter (marker-specific)
    LENGTH_FILTER(
        CHIMERA_REMOVAL.out.nochimera,
        params.marker,
        markerParams.min_length,
        markerParams.max_length
    )
    
    // Step 8: Generate ASV table
    GENERATE_ASV_TABLE(
        LENGTH_FILTER.out.asvs,
        CONCATENATE.out.concatenated,
        params.marker
    )
    
    // Step 9: BLAST taxonomic assignment
    if (!params.skip_blast) {
        BLAST_TAXONOMY(
            LENGTH_FILTER.out.asvs,
            params.marker
        )
        
        // Step 10: SSS filtering (optional)
        if (!params.skip_sss_filter && params.sss.apply_filter) {
            FILTER_SSS(
                BLAST_TAXONOMY.out.blast,
                GENERATE_ASV_TABLE.out.table,
                params.marker,
                params.sss.min_threshold
            )
        }
    }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    ========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output dir: ${params.outdir}
    ========================================
    """
}
