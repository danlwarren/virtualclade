#' Takes a vc.species object, a raster template, and a buffer width.  Breaks combined range apart into N clusters using kmeans clustering, and then assigns each cluster to the species whose PA centroid is closest to that cluster. Then it buffers the intersection of the cluster and the species PA out by buffer.width, and multiplies by the suitability raster.  That gives us the suitability of habitat within the actual range of the species.  Everything is then stuffed into a vc.clade object.
#'
#' @param x A vc.clade object to allopatrify
#' @param buffer.width Controls the post-allpatrification expansion of species' ranges.  Higher numbers create more range overlap.
#' @param plot Controls whether to print plots when generating PA data
#' @param split.cols Vector of column names to use for partitioning species occurrences.
#' @param beta Beta parameter past to convertToPA from virtual species
#' @param alpha Alpha parameter past to convertToPA from virtual species
#' @param env Optional stack of environmental layers, which will eventually be used to partition species based on environment as well as geography.
#' @param npres Number of presence points per species to use for partitioning ranges.
#' @param min.suitable Minimum number of "suitable" grid cells in a species range
#' @param nreps Number of replicates to try for permuting ranges to match phylogeny
#'
#' @export allopatrify
#'
#' @return A vc.clade object with ranges added to the "actual range" variable for each species

allopatrify <- function(x, buffer.width = 1, plot=TRUE, split.cols = c("lon", "lat"),
                        beta = 0.5, alpha = -0.00007, env = NA, npres = NA, min.suitable = 10,
                        nreps = 1000){

   # Get the species name template from the clade object
   species.var <- gsub(".1", "", names(x$species)[1])

   print("Getting pa data...")

   pa <-  clade.pa(x, sample.source = "suitab.raster",
                   method = "probability", beta = beta,
                   alpha = alpha, plot = plot, npres = npres)

   # Start building the output df
   pa.table <-  data.frame(pa$pa.table)

   # Set up raster stacks
   raster.template <- x$species[[1]]$virtualspecies$suitab.raster
   raster.template[!is.na(raster.template)] <- 0
   cluster.stack <- stack(raster.template)
   suit.stack <- stack(raster.template)
   pa.stack <- stack(raster.template)
   sample.stack <- stack(raster.template)

   # Get count of presence cells for each species, sorting smallest to largest
   breadth <- sort(table(pa.table[species.var]))

   if(inherits(env, c("raster", "RasterBrick", "RasterStack"))){
      env.cols <- names(env)[names(env) %in% split.cols]
      if(length(env.cols) > 0){
         pa.table <- cbind(pa.table,
                           raster::extract(env[[env.cols]], pa.table[,1:2]))
      }
   }

   message("allopatrifying...")

   # Get kmeans clusters
   unique.rows <- dplyr::distinct(pa.table, across(all_of(split.cols)), .keep_all = TRUE)
   clusters <- kmeans(unique.rows[,split.cols], length(unique(pa.table[,species.var])), nstart=1, iter.max=1000)
   unique.clusters <- cbind(unique.rows, clusters$cluster)
   colnames(unique.clusters)[ncol(unique.clusters)] <- "cluster"

   # Get centroids of clusters
   centroids <- plyr::ddply(unique.clusters[,c("cluster", split.cols)], "cluster", plyr::numcolwise(mean))

   # This is just for a quick prototype, might want to make it more general eventually
   cluster.raster <- rasterize(unique.clusters[,c("lon", "lat")], raster.template, field=unique.clusters[,"cluster"])
   raster::plot(cluster.raster, main = "Starting clusters")

   pa.table$cluster <- extract(cluster.raster, pa.table[,c("lon", "lat")])

   centroid.dist <- as.matrix(dist(centroids[,2:3], diag = TRUE, upper = TRUE))
   phylo.dist <- cophenetic(x$tree)

   suitable.cells <- pa.table %>%
      group_by(species, cluster) %>%
      summarise(n = n()) %>%
      as.data.frame()

   best.perm <- centroid.dist
   max.cor <- 0
   counter <- 1

   pb <- progress::progress_bar$new(total = nreps)

   while(counter < nreps){

      # Getting a possible random order
      inds <- sample(1:ncol(centroid.dist))

      # Creating an empty vector to store # suitable
      n.suitable <- numeric(length(inds))

      # Counting how many suitable (presence) grid cells each species has under this solution
      for(j in 1:length(inds)){
         this.species <- rownames(phylo.dist)[j]
         this.cluster <- rownames(centroid.dist)[inds[j]]
         this.suitable <- suitable.cells %>%
            filter(species == this.species, cluster == this.cluster) %>%
            select(n)
         n.suitable[j] <- as.numeric(this.suitable)
      }

      # Rejecting if any species doesn't have enough cells
      if(any(n.suitable < min.suitable) | any(is.na(n.suitable))){
         next
      }

      # Progress bar only ticks if we have enough suitable cells for each species
      pb$tick()

      # At this point we know the solution has enough cells per species, now we are
      # permuting the centroid distance matrix to see how well it matches the tree
      this.perm <- centroid.dist[inds,inds]
      this.cor <- cor(as.numeric(this.perm), as.numeric(phylo.dist), method = "spearman")

      # If this is the best solution so far we store it and update the max.cor value
      if(this.cor > max.cor){
         max.cor <- this.cor
         best.perm <- this.perm
      }
      counter <- counter + 1
   }

   cluster.assignments <- data.frame(species = rownames(phylo.dist),
                                     cluster = as.numeric(rownames(best.perm)))




   ###### THIS IS WHERE I STOPPED FOR A BIT!!!

   # Making a unique raster for each cluster, putting them in a stack
   for(i in cluster.assignments$cluster){
      cluster.stack <- addLayer(cluster.stack, cluster.raster == i)
   }
   cluster.stack <- dropLayer(cluster.stack, 1)
   names(cluster.stack) <- cluster.assignments$species

   # Making a raster of suitable habitat for each species, putting it in a stack
   for(i in cluster.assignments$species){
      suit.stack <- addLayer(suit.stack,
                             rasterize(pa.table[pa.table[species.var] == i,c("lon", "lat")], suit.stack, field = 1))
   }
   suit.stack <- dropLayer(suit.stack, 1)
   names(suit.stack) <- cluster.assignments$species

   if(plot == TRUE){
      raster::plot(pa.stack)
   }

   message("\n\nsympatrifying...")

   # Here's where we're generating sympatry.  In the original sim I had this set up
   # to produce a 1 everywhere that was present in the PA raster, the true suitability
   # in any grid cell that was within the buffer but NOT present in the PA raster,
   # and zero everywhere else.  I can't recall why I did it that way though. I think
   # for simulation purposes it might be better to have the suitability at every grid cell
   # in the range and 0 elsewhere

   for(i in cluster.assignments$species){
      temp.pa <- cluster.stack[[i]]
      temp.pa[temp.pa < 1 ] <- NA
      buffered.suitability <- raster::buffer(temp.pa, buffer.width) *
         x$species[[i]]$virtualspecies$suitab.raster
      buffered.suitability[is.na(buffered.suitability)] <- 0
      this.range <- raster.template +
         buffered.suitability
      this.range[this.range > 1] <- 1

      sample.stack[[i]] <- this.range
   }
   sample.stack <- dropLayer(sample.stack, 1)

   if(plot == TRUE){
      raster::plot(sample.stack)
   }

   # Creating a vc.clade object for output
   output <- x

   for(i in names(sample.stack)){
      output$species[[i]]$actual.range <- sample.stack[[i]]
   }

   # Might want to put more stuff into a list and return that
   return(output)
}
