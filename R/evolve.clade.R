#' Evolves a clade of virtual species from a single virtualspecies object.
#'
#' @param root.species A virtualspecies object to be used as the common ancestor of all simulated species.  Note that this currently has to be one that was generated with approach = "response".
#' @param ntaxa Number of tip taxa desired.
#' @param env A stack of environmental rasters.
#' @param rate The Brownian rate parameter controlling the rate of niche evolution.  Since different niche axes will have different parameter values, this rate is expressed as as a proportion of the value of each parameter at the root node. So if rate = 0.5 and the value of a parameter for the root species is 4, the sigma parameter for ape's rTraitCont simulation will be 0.5 * 4 = 2.
#' @param min.sd The minimum standard deviation for each parameter, as a proportion of the parameter value at the root.  Essentially this sets the minimum possible niche width.
#' @param tip.prefix The prefix to use for naming species.
#'
#' @examples
#' my.root <- generateRandomSp(euro.worldclim[[c(1,12)]])
#' my.clade <- evolve.clade(my.root, 10, euro.worldclim)
#'
#' @export evolve.clade
#'
#' @return A list containing the root species, the tip species, and a phylogeny.  Also contains the arguments used to run the simulation.


evolve.clade <- function(root.species, ntaxa, env, rate = 0.2, min.sd = .5, tip.prefix = "species", max.tries = 10, ...){

   if("pca" %in% names(root.species$details)){
      stop("evolve.clade only works with virtual species generated using approach = 'response'")
   }

   # Grabbing the params from our root species to evolve them
   root.species <- unpack.species(root.species)
   root.params <- root.species$details$parameters

   # Creating names for our tips
   tip.names <- paste(tip.prefix, 1:ntaxa, sep=".")

   # Make a random tree
   tree <- ape::rtree(ntaxa, tip.label = tip.names)

   # Get the names of the env vars, trim the stack to just those (virtualspecies requirement)
   env.vars <- root.species$details$variables
   env <- env[[env.vars]]

   param.list <- list()

   # Creating a copy of root.params for each species, just to get the structure.
   # we'll modify these from the sims, and then we'll use them to build new species
   # from generateSpFromFun in virtualspecies
   for(i in tip.names){
      param.list[[i]] <- root.params
   }

   tries = 1

   while(tries <= max.tries){
      for(i in names(root.params)){
         # i is iterating over predictor variables
         for(j in names(root.params[[i]]$args)){

            if(grepl("sd", j)){
               # Brownian, but constraining the sd to be > a proportion of the root value
               this.root.value <- root.params[[i]][["args"]][j]


               # Not entirely clear why I have to make this.root.value numeric explicitly, but
               # apparently I do
               this.trait <- ape::rTraitCont(tree,
                                             root.value = this.root.value,
                                             theta = this.root.value,
                                             sigma = rate * as.numeric(this.root.value))

               for(k in tip.names){
                  param.list[[k]][[i]][["args"]][j] <- max(abs(this.trait[k]), min.sd * this.root.value)
               }
            }
            else{

               this.root.value <- root.params[[i]][["args"]][j]


               # Not entirely clear why I have to make this.root.value numeric explicitly, but
               # apparently I do
               this.trait <- ape::rTraitCont(tree,
                                             root.value = this.root.value,
                                             theta = this.root.value,
                                             sigma = rate * as.numeric(this.root.value))

               for(k in tip.names){
                  param.list[[k]][[i]][["args"]][j] <- this.trait[k]
               }
            }
         }
      }

      species <- list()
      for(i in names(param.list)){
         species[[i]] <- vc.species(virtualspecies =  suppressMessages(generateSpFromFun(raster.stack = env,
                                                                        parameters = param.list[[i]],
                                                                        plot = TRUE)),
                                    species.name = i)
      }

      # Check to see if any species have gone extinct
      if(any(is.na(sapply(species, function(x) maxValue(x$virtualspecies$suitab.raster))))){
         tries <- tries + 1
         message(paste0("Some species went extinct, attempt ", tries))
      } else {
         break
      }
   }


   root.species <- vc.species(virtualspecies = root.species, species.name = "root.species")


   output <- vc.clade(species = species,
                         tree = tree,
                         root.species = root.species,
                         sim.args = as.list(match.call()))

   class(output) <- c("list", "vc.clade")

   return(output)
}







