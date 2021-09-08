#' Plots parameter values for a set of species in a virtualclade object.  Returns both histograms of the distribution of parameter values and a set of trees with parameter values mapped on them.
#'
#' @param clade A virtualclade object
#'
#' @examples
#' my.root <- generateRandomSp(euro.worldclim[[c(1,12)]])
#' my.clade <- evolve.clade(my.root, 10, euro.worldclim)
#' parameter.plots(my.clade)
#'
#' @export parameter.plots
#'
#' @return A list containing the plotting data frame and a series of plots.


parameter.plots <- function(clade){

  # Get a vector of parameter names
  root.species <- unpack.species(clade$root.species$virtualspecies)
  root.params <- root.species$details$parameters

  param.df <- data.frame(species = character(0),
                         parameter = character(0),
                         value = numeric(0))



  for(i in names(root.params)){
    # i is iterating over predictor variables
    for(j in names(root.params[[i]]$args)){
      root.value <- root.params[[i]]$args[[j]]
      this.param <- paste0(i, ".", j)
      param.df <- rbind(param.df, c("root", this.param, root.value))

      for(k in names(clade$species)){
        this.value <- clade$species[[k]]$virtualspecies$details$parameters[[i]]$args[[j]]
        param.df <- rbind(param.df, c(k, this.param, this.value))
      }
    }
  }

  colnames(param.df) <- c("species", "parameter", "value")
  param.df$value <- as.numeric(param.df$value)

  param.hist <- qplot(value, data = param.df, fill = parameter, geom = "histogram") +
    stat_bin(bins = length(clade$species)) +
    facet_wrap(. ~ parameter, scales = "free") +
    theme_bw()

  treeplots <- list()

  for(i in unique(param.df$parameter)){
    trait <- param.df[which(param.df$parameter == i),]
    trait <- trait[trait$species %in% clade$tree$tip.label,]
    trait.vec <- as.numeric(trait$value)
    names(trait.vec) <- trait$species

    obj <- phytools::contMap(clade$tree, trait.vec, plot=FALSE)
    obj <- phytools::setMap(obj, invert=TRUE)
    plot(obj, outline=FALSE, lwd=c(3,7), leg.txt=i)
    treeplots[[i]] <- recordPlot()
  }





  output <- list(param.hist = param.hist,
                 treeplots = treeplots,
                 param.df = param.df)
  return(output)
}

