nextflow.enable.dsl = 2

include { initOptions } from '../functions'

process DL_VIRAL_PROTEIN_DB {
    publishDir params.dbdir, mode: 'copy'
    container = 'nakor/virus_extraction'

    output:
    path 'refseq_viral_proteins.dmnd'

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.protein.faa.gz
    cat viral.*.protein.faa.gz | unpigz -p $task.cpus > viral.protein.faa
    diamond makedb --threads $task.cpus --in viral.protein.faa --db refseq_viral_proteins
    """    
    
}

process DIAMOND_BLASTX {
    tag {"${meta.id}"}
    publishDir params.outdir+"/diamond", mode: "copy"
    container = 'nakor/virus_extraction'
    
    input:
    tuple val(meta), path(fasta)
    path db
    val options

    output:
    tuple val(meta), path('*'), emit: 'all'
    tuple val(meta), path('diamond_blastx_contigs*.txt'), emit: 'ctg_ids'    
    
    script:
    def ioptions = initOptions(options)
    """
    diamond blastx $ioptions.args -d ${db} -q ${fasta} -o diamond_matches_on_refseq_${meta.id}.m8 --threads $task.cpus
    cut -f1 diamond_matches_on_refseq_${meta.id}.m8 | sort | uniq -c | sort -rnk1 \
        | awk '{OFS=","}{print \$2,\$1}' > diamond_ctg_hit_counts_${meta.id}.txt
    cut -d, -f1 diamond_ctg_hit_counts_${meta.id}.txt > diamond_blastx_contigs_${meta.id}.txt
    """
}

workflow diamond_blastx {
    take:
    contigs // tuple (meta, fasta)
    options

    main:
    vir_prot_db = file(params.vir_prot_db)
    if(!vir_prot_db.exists()) {
        vir_prot_db = DL_VIRAL_PROTEIN_DB()
    }
    DIAMOND_BLASTX(contigs, vir_prot_db, options)

    emit:
    all = DIAMOND_BLASTX.out.all
    ctg_ids = DIAMOND_BLASTX.out.ctg_ids
}

workflow test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    diamond_blastx(fasta, [:])
}
