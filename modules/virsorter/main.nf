nextflow.enable.dsl = 2

include { initOptions } from '../functions'

process DL_VIRSORTER_DB {
    publishDir params.dbdir, mode: 'copy'

    output:
    path 'vs_db'

    script:
    """
    wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
    tar -xvzf virsorter-data-v2.tar.gz
    mv virsorter-data vs_db
    """    
}

process VIRSORTER {
    tag {"${meta.id}"}
    container 'nakor/virsorter'
    
    publishDir params.outdir+"/virsorter", mode: "copy"

    input:
    tuple val(meta), path(fasta)
    path vs_db
    val options

    output:
    tuple val(meta), path('vs_out'), emit: all
    tuple val(meta), path ('virsorter_contigs*.txt'), emit: ctg_ids
    
    script:
    def ioptions = initOptions(options)
    """
    wrapper_phage_contigs_sorter_iPlant.pl $ioptions.args \
        -f $fasta \
        --ncpu $task.cpus \
        --wdir vs_out \
        --data-dir $vs_db --db 1
    
    grep "^VIRSorter" vs_out/VIRSorter_global-phage-signal.csv \
        | cut -d, -f1 \
        | sed 's/VIRSorter_//' \
        | sed 's/-circular//' \
        > virsorter_contigs_${meta.id}.txt
    """
}

process DL_VIRSORTER2_DB {
    publishDir params.dbdir, mode: 'copy'
    container 'nakor/virsorter'

    output:
    path 'virsorter2_db'

    script:
    """
    virsorter setup -d virsorter2_db -j $task.cpus
    """    
}

process VIRSORTER2 {
    tag {"${meta.id}"}
    publishDir params.outdir+"/virsorter2", mode: "copy"
    publishDir params.outdir, mode: "copy", pattern: "virsorter2_contigs*.txt"
    container 'nakor/virsorter'

    input:
    tuple val(meta), path(fasta)
    path vs2_db
    val options

    output:
    tuple val(meta), path('final-viral-score.tsv'), emit: all
    tuple val(meta), path ('virsorter2_contigs*.txt'), emit: ctg_ids
    
    script:
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env bash

    virsorter run $ioptions.args  -w out -i $fasta -j $task.cpus -d $vs2_db
    mv out/final-viral-score.tsv .
    tail -n+2 final-viral-score.tsv | cut -f1 | cut -d'|' -f1 > virsorter2_contigs_${meta.id}.txt
    """
}


workflow virsorter {
    take:
    contigs // tuple (meta, fasta)
    options

    main:
    vs_db = file(params.vs_db)
    if(!vs_db.exists()) {
        vs_db = DL_VIRSORTER_DB()
    }
    VIRSORTER(contigs, vs_db, options)

    emit:
    all = VIRSORTER.out.all
    ctg_ids = VIRSORTER.out.ctg_ids
}

workflow virsorter2 {
    take:
    contigs // tuple (meta, fasta)
    options

    main:
    vs2_db = file(params.vs2_db)
    if(!vs2_db.exists()) {
        vs2_db = DL_VIRSORTER2_DB()
    }
    VIRSORTER2(contigs, vs2_db, options)

    emit:
    all = VIRSORTER2.out.all
    ctg_ids = VIRSORTER2.out.ctg_ids
}

workflow vs_test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    virsorter(fasta, [:])
}

workflow vs2_test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    virsorter2(fasta, [:]) // doesn't work yet
}

workflow test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    virsorter(fasta, [:])
    // virsorter2(fasta, [:]) // doesn't work yet
}
