#' Defining a class for vc.clade, which is a container for clades of simulated virtual species and related items.
#'
#' @param species A list of vc.species objects.
#' @param tree A tree showing the relationships between the species.
#' @param root.species A vc.species object for the root species.
#' @param sim.args The arguments used to call evolve.clade for the initial simulation.
#'
#' @export vc.clade
#'
#' @return A vc.clade object.


vc.clade <- function(species = NA, tree = NA, root.species = NA, sim.args = NA){

   # Checking classes of input args.  The isTRUE stuff is needed because R doesn't
   # know how to do is.na on raster data, so it was barfing and error when a raster
   # was passed in.

   if(!isTRUE(is.na(species))){

      # Checking to see if species is a list
      if(!"list" %in% class(species)){
         print("Argument species requires a list of vc.species objects")
      }

      # This if statement is asking whether any of the list elements don't have
      # vc.species in their class definitions
      if(any(unlist(lapply(species, function(x) !"vc.species" %in% class(x))))){
         print("The following objects in the species list do not appear to be vc.species objects:")
         print(names(which(unlist(lapply(species, function(x) !"vc.species" %in% class(x))))))
      }

   }

   if(!isTRUE(is.na(tree))){
      # Checking to see if species is a list
      if(!"phylo" %in% class(tree)){
         print("Argument tree requires a phylo object")
      }
   }

   if(!isTRUE(is.na(root.species))){
      # Checking to see if species is a list
      if(!"vc.species" %in% class(root.species)){
         print("Argument root.species requires a vc.species object")
      }
   }

   output <- list(species = species,
                  tree = tree,
                  root.species = root.species,
                  sim.args = sim.args)

   class(output) <- c("list", "vc.clade")

   return(output)
}

plot.vc.clade <- function(this.clade){
   lapply
}
