
include { TBL_2_JSON as TBL_2_JSON_1V1      } from '../../modules/local/tbl_2_json.nf'
include { TBL_2_JSON as TBL_2_JSON_1VM      } from '../../modules/local/tbl_2_json.nf'
include { COMBINE_JSON                      } from '../../modules/local/combine_json.nf'
include { VOTE_BEST_MATCH                   } from '../../modules/local/vote_bm.nf'

workflow CONSENSUS_VOTE {

    take:
    ch_ortho_l     // channel: [meta, file]  — Orthologer  (1v1, per reference)
    ch_ortho_f     // channel: [meta, file]  — OrthoFinder (1v1)
    ch_eggnog      // value:   file          — EggNOG-mapper (1vM)
    ch_tree        // value:   file          — TreeGrafter (1vM)
    ch_synteny     // channel: [meta, file]  — MCScan synteny (1v1, per reference; empty if run_synteny=false)
    ch_query_seqs  // value:   path          — query protein FASTA (from TransDecoder)

    main:

    ch_1v1_tools  = Channel.empty()
    ch_1vm_tools  = Channel.empty()
    ch_final_jsons = Channel.empty()

    // --- 1v1 evidence sources ---

    ch_ortho_l.map  { meta, file -> [ meta + [ tool: "ortho_l"  ], file ] }
        .set { ch_ortho_l_comb }
    ch_1v1_tools.mix(ch_ortho_l_comb)
        .set { ch_1v1_tools }

    ch_ortho_f.map  { meta, file -> [ meta + [ tool: "ortho_f"  ], file ] }
        .set { ch_ortho_f_comb }
    ch_1v1_tools.mix(ch_ortho_f_comb)
        .set { ch_1v1_tools }

    // Synteny is 1v1 (one best syntenic partner per query gene per reference species)
    ch_synteny.map  { meta, file -> [ meta + [ tool: "synteny"  ], file ] }
        .set { ch_synteny_comb }
    ch_1v1_tools.mix(ch_synteny_comb)
        .set { ch_1v1_tools }

    // --- 1vM evidence sources ---

    ch_1vm_tools.mix( ch_eggnog.map { file -> [ [ tool: "eggnog" ], file ] } )
        .set { ch_1vm_tools }
    ch_1vm_tools.mix( ch_tree.map   { file -> [ [ tool: "tree"   ], file ] } )
        .set { ch_1vm_tools }

    // --- Convert to JSON ---

    TBL_2_JSON_1V1(ch_1v1_tools)

    TBL_2_JSON_1VM(ch_1vm_tools)
    ch_final_jsons.mix(TBL_2_JSON_1VM.out.map { meta, file -> file })
        .set { ch_final_jsons }

    ch_all_jsons = TBL_2_JSON_1V1.out.json
    ch_all_jsons.map { meta, json -> [ meta.tool, json ] }
        .groupTuple()
        .map { tool, jsons -> [ [ tool: tool ], jsons ] }
        .set { ch_all_jsons }

    COMBINE_JSON(ch_all_jsons)
    ch_final_jsons.mix(COMBINE_JSON.out.map { meta, file -> file })
        .set { ch_final_jsons }

    ch_final_jsons.collect().set { ch_final_jsons }

    VOTE_BEST_MATCH(
        ch_final_jsons,
        ch_query_seqs,
    )

    emit:
    orthologs = VOTE_BEST_MATCH.out.combined_orthologs

}
