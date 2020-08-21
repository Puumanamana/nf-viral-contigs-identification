nextflow.enable.dsl = 2

include {viralverify} from './modules/viralverify/main.nf' addParams(outdir: "${params.outdir}/virus_detection")
include {virsorter; virsorter2} from './modules/virsorter/main.nf' addParams(outdir: "${params.outdir}/virus_detection")
include {virfinder} from './modules/virfinder/main.nf' addParams(outdir: "${params.outdir}/virus_detection")
include {vibrant} from './modules/vibrant/main.nf' addParams(outdir: "${params.outdir}/virus_detection")
include {diamond_blastx} from './modules/diamond/main.nf' addParams(outdir: "${params.outdir}/virus_detection")
include {CONTIGS_FROM_IDS; SUMMARIZE_CALLS} from './modules/misc/main.nf' addParams(outdir: "${params.outdir}/virus_detection")


def init_all(Map args) {
    for (key in ['virsorter', 'virfinder', 'vibrant', 'viralverify', 'diamond']) {
        if (!args.containsKey("${key}_args")) {
            args["${key}_args"] = ""
        }
    }
}

workflow {
    assembly = Channel.fromPath(params.assembly).map{[[id: it.getSimpleName()], it]}

    init_all(params)
    
    // run multiple softwares
    virsorter(assembly, [args: params.virsorter_args])
    virfinder(assembly, [args: params.virfinder_args])
    viralverify(assembly, [args: params.viralverify_args])
    vibrant(assembly, [args: params.vibrant_args])
    diamond_blastx(assembly, [args: params.diamond_args])

    virsorter.out.ctg_ids.mix(
        virfinder.out.ctg_ids,
        viralverify.out.ctg_ids,
        vibrant.out.ctg_ids,
        diamond_blastx.out.ctg_ids
    ).set{all_ids}

    SUMMARIZE_CALLS(
        [samples: assembly.collect{it[0].id},
         tools: ['virsorter', 'virfinder', 'viralverify', 'vibrant', 'diamond_blastx']],
        all_ids.collect{it[1]}
    )
}
