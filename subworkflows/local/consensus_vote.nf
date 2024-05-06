
include { TBL_2_JSON as TBL_2_JSON_1V1      } from '../../modules/local/tbl_2_json.nf'
include { TBL_2_JSON as TBL_2_JSON_1VM      } from '../../modules/local/tbl_2_json.nf'
include { COMBINE_JSON                      } from '../../modules/local/combine_json.nf'
include { VOTE_BEST_MATCH                   } from '../../modules/local/vote_bm.nf'

workflow CONSENSUS_VOTE {

     take:
     ch_ortho_l     // tuple: [meta, orthologer outputs]
     ch_ortho_f     // tuple: [meta, orthofinder outputs]
     ch_eggnog      // value:    eggnogmapper output
     ch_tree        // value:    treegrafter output

     main:
     
     ch_1v1_tools = Channel.empty()
     ch_1vm_tools = Channel.empty()

     ch_ortho_l.map { meta, file -> [ meta + [ tool: "ortho_l" ], file ] }
     .set { ch_ortho_l_comb }
     ch_1v1_tools.mix(ch_ortho_l_comb)
     .set { ch_1v1_tools }
     
     ch_ortho_f.map { meta, file -> [ meta + [ tool: "ortho_f" ], file ] }  
     .set { ch_ortho_f_comb }
     ch_1v1_tools.mix(ch_ortho_f_comb)
     .set { ch_1v1_tools }
       
     ch_1vm_tools.mix( ch_eggnog.map { file -> [ [ tool: "eggnog" ], file ] } ).set { ch_1vm_tools }
     ch_1vm_tools.mix( ch_tree.map { file -> [ [ tool: "tree" ], file ] } ).set { ch_1vm_tools }
     
     
     TBL_2_JSON_1V1(ch_1v1_tools)
     
     TBL_2_JSON_1VM(ch_1vm_tools)
     
     ch_all_jsons = TBL_2_JSON_1V1.out.json
     ch_all_jsons.map { meta, json -> [ meta.tool, json ] }
     .groupTuple()
     .map { meta, jsons -> [ [ tool: meta ], jsons ] }
     .set { ch_all_jsons }

     ch_all_jsons.view()     
     COMBINE_JSON(ch_all_jsons)
     
     //VOTE_BEST_MATCH(
             //COMBINE_JSON.out.ortho_f,
            //COMBINE_JSON.out.ortho_l,
           //TBL_2_JSON_1VM.out.eggnog,
          // TBL_2_JSON_1VM.out.tree
         // )
     
    
     

     emit:
     tbl_test = TBL_2_JSON_1V1.out   



}

