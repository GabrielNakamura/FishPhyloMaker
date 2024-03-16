
#' Calculate PD deficit at species and insertion levels
#'
#' @param phylo a phylo object
#' @param data data frame with species and insertion levels. This is the output of FishPhyloMaker function
#'
#' @return a list with two tibbles. One with results of deficit for all species and another with summary results for each category of insertion
#' @importFrom magrittr "%>%"
#' @return a list containing two data frames. 
#' 
#'     - a long data frame (pd_deficit_long) containing all values of deficit and total phylogenetic diversity
#'         accounted by each terminal branch of the tree.
#'     
#'     - a summary data frame with values of total PD in each category of insertion (including for those species present in tree)
#'         and a proportion of each category for the total terminal PD.
#' @export
#'
#' @examples
#' 
PD_deficit_2 <- 
  function(phylo, data){
    tree_tib <- tidytree::as_tibble(phylo) # change to a tibble object  
    
    # join with insertion data
    tree_tib2 <- 
      tree_tib %>% 
      right_join(data, by = c("label" = "s")) 
    
    # Summary stats for whole table
    tibble_res_all <- 
      tree_tib2 %>% 
      group_by(insertions) %>% 
      add_count(name = "n.insertion.category") %>% 
      mutate(pd.total.category = sum(branch.length), 
             pd.relative.category = sum(branch.length)/sum(phylo$edge.length), 
             pd.relative.spp = branch.length/sum(phylo$edge.length)) %>% 
      rename(species = "label", 
             order = "o",
             family = "f") %>% 
      select(species, order, family, insertions, branch.length, n.insertion.category, pd.total.category, pd.relative.category, pd.relative.spp)
    
    # preparing the output with summary stats
    summary_table <- 
      tibble_res_all %>% 
      distinct(insertions, n.insertion.category, pd.total.category, pd.relative.category)
    
    list_res <- vector(mode = "list")
    list_res$pd_deficit_long <- tibble_res_all
    list_res$summary_deficit <- summary_table
    return(list_res)
  }