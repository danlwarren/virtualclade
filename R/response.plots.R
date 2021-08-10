#' Plots responses to all environmental variables for a set of species in a virtualclade object
#'
#' @param clade A virtualclade object
#'
#' @examples
#' my.root <- generateRandomSp(euro.worldclim[[c(1,12)]])
#' my.clade <- evolve.clade(my.root, 10, euro.worldclim)
#' response.plots(my.clade)
#'
#' @export response.plots
#'
#' @return A list containing the plotting data frame and a ggplot object showing the response of each species to the environmental gradient.


response.plots <- function(clade){

  # Get a vector of variable names
  vars <- clade$root.species$virtualspecies$details$variables

  # Create an empty plot data frame.  We will add to this variable by variable
  # and species by species
  plot.df <- data.frame(species = character(0),
                        variable = character(0),
                        value = numeric(0),
                        suitability = numeric(0))



  for(i in vars){

    # First doing the root species
    parameters <- clade$root.species$virtualspecies$details$parameters

    # I'm borrowing the next little chunk of this code directly from virtualspecies
    var.seq <- seq(parameters[[i]]$min,
                   parameters[[i]]$max,
                   length = 100)
    this.pred <- do.call(match.fun(parameters[[i]]$fun), args = c(list(var.seq), parameters[[i]]$args))

    # Attaching the root species to the plot data frame
    plot.df <- rbind(plot.df, data.frame(species = rep("root.species", 100),
                                         variable = rep(i, 100),
                                         value = var.seq,
                                         suitability = this.pred))

    # Now we're going to step through the descendent species and do the same for them
    for(j in names(clade$species)){

      # Get the params for this species and variable
      parameters <- clade$species[[j]]$virtualspecies$details$parameters

      var.seq <- seq(parameters[[i]]$min,
                     parameters[[i]]$max,
                     length = 100)
      this.pred <- do.call(match.fun(parameters[[i]]$fun), args = c(list(var.seq), parameters[[i]]$args))

      # Attach to plot data frame
      plot.df <- rbind(plot.df, data.frame(species = rep(j, 100),
                                           variable = rep(i, 100),
                                           value = var.seq,
                                           suitability = this.pred))
    }
  }

  # Scale the suitability values so they're (0,1) across each variable
  plot.df <- plot.df %>%
    group_by(variable) %>%
    mutate(suitability = (suitability - min(suitability)) / (max(suitability) - min(suitability)))

  response.plot <- qplot(value, suitability, data = plot.df,
                         color = species, geom = "line") +
    facet_grid(~variable, scales = "free") + theme_bw()

  output <- list(response.plot = response.plot,
                 plot.df = plot.df)
  return(output)
}

