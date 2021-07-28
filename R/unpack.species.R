#' Takes a virtualspecies object and unlists all of the args in the niche parameters.
#' If you don't do this after generating random species, ape will barf and die horribly.
#'
#' @param species A virtualspecies object
#'
#'
#' @return A virtualspecies object with its guts unpacked slightly

unpack.species <- function (species){
  for(i in species$details$variables){
    species$details$parameters[[i]]$args <- unlist(species$details$parameters[[i]]$args)
  }
  return(species)
}
