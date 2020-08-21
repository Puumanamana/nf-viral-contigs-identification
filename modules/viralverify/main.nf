nextflow.enable.dsl = 2

include { initOptions } from '../functions'

process DL_PFAM_DB {
    tag {"download_pfam_db"}
    publishDir params.dbdir, mode: 'copy'
    container 'nakor/virus_extraction'
    
    output:
    file('Pfam-A.hmm')

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz 
    unpigz -p $task.cpus Pfam-A.hmm.gz
    """
}

process VIRALVERIFY {
    tag {"${meta.id}"}
    publishDir params.outdir+"/viralverify", mode: "copy"
    container 'nakor/virus_extraction'
    
	input:
    tuple val(meta), path(fasta)
    path pfam_db
    val options

	output:
    tuple val(meta), path('Prediction_results_fasta/*'), emit: all
    tuple val(meta), path('viralverify_contigs*.txt'), emit: ctg_ids

    script:
    def ioptions = initOptions(options)
    """
    viralverify.py $ioptions.args -f ${fasta} -o ./ --hmm ${pfam_db} -t $task.cpus
    cat Prediction_results_fasta/*virus* | grep '^>' | cut -c2- > viralverify_contigs_${meta.id}.txt
    """
}


workflow viralverify {
    take:
    contigs
    options

    main:
    pfam_db = file(params.pfam_db)
    if(!pfam_db.exists()) {
        pfam_db = DL_PFAM_DB()
    }
    VIRALVERIFY(contigs, pfam_db, options)

    emit:
    all = VIRALVERIFY.out.all
    ctg_ids = VIRALVERIFY.out.ctg_ids
}

workflow test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    viralverify(fasta, [:])
}
