nextflow.enable.dsl = 2

include { initOptions } from '../functions'

process DL_VIBRANT_DB {
    tag {"download_vibrant_db"}
    publishDir params.dbdir, mode: 'copy'
    container 'nakor/vibrant'
    
    output:
    path "$db_name/databases", emit: db
    path "$db_name/files", emit: files

    script:
    def db_name = params.vibrant_db.tokenize('/')[-1]
    """
    mkdir -p "${db_name}/databases" && cd "${db_name}/databases"
    vibrant_dir=\$(dirname \$(which VIBRANT_setup.py))
    cp -r "\$vibrant_dir"/profile_names .
    cp -r "\$vibrant_dir"/../files ..

    wget -qO- http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz | tar xz
    wget -q ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz \
        && gunzip Pfam-A.hmm.gz
    wget -qO- ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz | tar xz

    for v in VOG*.hmm; do 
	    cat \$v >> vog_temp.HMM
    done

    for k in profiles/K*.hmm; do 
	    cat \$k >> kegg_temp.HMM
    done

    rm -r VOG0*.hmm VOG1*.hmm VOG2*.hmm profiles

    hmmfetch -o VOGDB94_phage.HMM -f vog_temp.HMM profile_names/VIBRANT_vog_profiles.txt >> VIBRANT_setup.log
    hmmfetch -o KEGG_profiles_prokaryotes.HMM -f kegg_temp.HMM profile_names/VIBRANT_kegg_profiles.txt >> VIBRANT_setup.log
    
    rm vog_temp.HMM kegg_temp.HMM
    mv Pfam-A.hmm Pfam-A_v32.HMM

    hmmpress VOGDB94_phage.HMM >> VIBRANT_setup.log
    hmmpress KEGG_profiles_prokaryotes.HMM >> VIBRANT_setup.log
    hmmpress Pfam-A_v32.HMM >> VIBRANT_setup.log
    """
}

process VIBRANT {
    tag {"${meta.id}"}
    publishDir params.outdir+"/vibrant", mode: "copy"
    container 'nakor/vibrant'
    
	input:
    tuple val(meta), path(fasta)
    path vibrant_db
    path vibrant_files
    val options

	output:
    tuple val(meta), path('*'), emit: all
    tuple val(meta), path('vibrant_contigs*.txt'), emit: ctg_ids

    script:
    def ioptions = initOptions(options)
    def radical = "${fasta.getSimpleName()}"
    def result_path = "VIBRANT_${radical}/VIBRANT_results_${radical}/VIBRANT_machine_${radical}.tsv"
    """
    VIBRANT_run.py $ioptions.args -i $fasta -t $task.cpus -d $vibrant_db -m $vibrant_files
    awk -F"\\t" '\$2=="virus"' $result_path | cut -d ' ' -f 1 >> vibrant_contigs_${meta.id}.txt
    """
}


workflow vibrant {
    take:
    contigs
    options

    main:
    vibrant_data = [db: file("${params.vibrant_db}/databases"), files: file("${params.vibrant_db}/files")]
    if(!vibrant_data.db.exists()) {
        vibrant_data = DL_VIBRANT_DB()
    }
    VIBRANT(contigs, vibrant_data.db, vibrant_data.files, options)

    emit:
    all = VIBRANT.out.all
    ctg_ids = VIBRANT.out.ctg_ids
}

workflow test {
    fasta = Channel.fromPath("$baseDir/../../test_data/sample*.fasta")
        .map{[[id: it.getSimpleName()], it]}
    vibrant(fasta, [:])
}
