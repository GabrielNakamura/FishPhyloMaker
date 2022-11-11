pb_names_family <- progress::progress_bar$new(format = "  Naming family nodes [:bar] :percent", 
                                              total = length(list_family), clear = FALSE, 
                                              width = 60, current = "<", incomplete = ">", 
                                              complete = ">")


if (progress.bar == TRUE) {
  pb_names_family_toadd <- progress::progress_bar$new(format = "  Naming monotipic family nodes [:bar] :percent", 
                                                      total = length(list_family_tobeaddnames), clear = FALSE, 
                                                      width = 60, current = "<", incomplete = ">", 
                                                      complete = ">")
}