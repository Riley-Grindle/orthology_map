



workflow {
 
    ch_test = Channel.of(0)


    ch_test.map { item -> "ZERO" }
    ch_test.view()

}
