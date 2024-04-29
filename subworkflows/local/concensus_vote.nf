\
include { TBL_2_JSON as TBL_2_JSON_1V1      } from '../../modules/local/tbl_2_json.nf'
include { TBL_2_JSON as TBL_2_JSON_1VM      } from '../../modules/local/tbl_2_json.nf'
include { COMBINE_JSON                      } from '../../modules/local/combine_json.nf'
include { VOTE_BEST_MATCH                   } from '../../modules/local/vote_bm.nf'

workflow CONCENSUS_VOTE {

     take:
     ch_ortho_l     // tuple: [meta, orthologer outputs]
     ch_ortho_f     // tuple: [meta, orthofinder outputs]
     ch_eggnog      // value:    eggnogmapper output
     ch_tree        // value:    treegrafter output

     main:
     
     ch_1v1_tools = Channel.empty()

     ch_ortho_l.map { meta, file -> [ meta + [ tool: "ortho_l" ], file ] }
     .set { ch_ortho_l_comb }
     ch_1v1_tools.mix(ch_ortho_l_comb)
     
     ch_ortho_f.map { meta, file -> [ meta + [ tool: "ortho_f" ], file ] }  
     .set { ch_ortho_f_comb }
     ch_1v1_tools.mix(ch_ortho_f_comb)
       
     ch_1vm_tools.mix( ch_eggnog.map { file -> [ [ tool: "eggnog" ], file ] } )
     ch_1vm_tools.mix( ch_tree.map { file -> [ [ tool: "tree" ], file ] } )

     TBL_2_JSON_1V1(ch_1v1_tools)
     
     TBL_2_JSON_1VM(ch_1vm_tools)
     
     
     //COMBINE_JSON(TBL_2_JSON_1V1.out.collect())
     
     //VOTE_BEST_MATCH(
           COMBINE_JSON.out.ortho_f,
           COMBINE_JSON.out.ortho_l,
           TBL_2_JSON_1VM.eggnog,
           TBL_2_JSON_1VM.tree
       )
     
    
     

     emit:
     peptide_fasta = ch_combined 
     versions      = ch_versions   



}

