#'  Takes a vc.species object, a raster template, and a buffer width.  Breaks combined range apart into N clusters using kmeans clustering, and then assigns each cluster to the species whose PA centroid is closest to that cluster. Then it buffers the intersection of the cluster and the species PA out by buffer.width, and multiplies by the suitability raster.  That gives us the suitability of habitat within the actual range of the species.  Everything is then stuffed into a vc.clade object.
#'
#'  @param x A vc.clade object to allopatrify
#'  @param buffer.width Controls the post-allpatrification expansion of species' ranges.  Higher numbers create more range overlap.
#'  @param plot Controls whether to print plots when generating PA data
#'  @param split.cols Vector of column names to use for partitioning species occurrences.
#'  @param beta Beta parameter past to convertToPA from virtual species
#'  @param alpha Alpha parameter past to convertToPA from virtual species
#'  @param env Optional stack of environmental layers, which will eventually be used to partition species based on environment as well as geography.
#'
#'  @export allopatrify
#'
#'  @return A vc.clade object with ranges added to the "actual range" variable for each species

allopatrify <- function(x, buffer.width = 1, plot=TRUE, split.cols = c("lon", "lat"), beta = 0.5, alpha = -0.00007, env = NA){

   # Get the species name template from the clade object
   species.var <- gsub(".1", "", names(x$species)[1])

   print("Getting pa data...")

   pa <-  clade.pa(x, sample.source = "suitab.raster", method = "probability", beta = beta, alpha = alpha, plot = plot)

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

   # Get centroids of suitable habitat
   #print(head(x.pa))
   centroids <- plyr::ddply(pa.table[,c(species.var, split.cols)], species.var, plyr::numcolwise(mean))

   # Get kmeans clusters
   unique.rows <- distinct(pa.table, across(all_of(split.cols)), .keep_all = TRUE)
   clusters <- kmeans(unique.rows[,split.cols], length(unique(pa.table[,species.var])), nstart=1, iter.max=1000)
   unique.clusters <- cbind(unique.rows, clusters$cluster)
   colnames(unique.clusters)[ncol(unique.clusters)] <- "cluster"

   # This is just for a quick prototype, might want to make it more general eventually
   cluster.raster <- rasterize(unique.clusters[,c("lon", "lat")], raster.template, field=unique.clusters[,"cluster"])
   raster::plot(cluster.raster)

   # Making a unique raster for each cluster, putting them in a stack
   for(i in unique(unique.clusters[,"cluster"])){
      cluster.stack <- addLayer(cluster.stack, cluster.raster == i)
   }
   cluster.stack <- dropLayer(cluster.stack, 1)
   names(cluster.stack) <- paste("species", unique(unique.clusters[,"cluster"]), sep = ".")

   # Making a raster of suitable habitat for each species, putting it in a stack
   for(i in unique(pa.table[,species.var])){
      suit.stack <- addLayer(suit.stack,
                             rasterize(pa.table[pa.table[species.var] == i,c("lon", "lat")], suit.stack, field = 1))
   }
   suit.stack <- dropLayer(suit.stack, 1)
   names(suit.stack) <- unique(pa.table[,species.var])

   print("allopatrifying...")

   # Stepping through from narrowest range to largest, putting species where they overlap with their cluster
   temp.cluster <- cluster.stack
   for(i in names(breadth)){
      #print(i)
      this.sp <- suit.stack[[i]]
      overlaps <- c()
      for(j in names(temp.cluster)){
         overlaps <- c(overlaps, cellStats(this.sp * temp.cluster[[j]], stat="sum"))
      }
      this.pa <- this.sp * temp.cluster[[which.max(overlaps)]]

      # Some species don't overlap with any remaining clusters, they get their original range back
      if(cellStats(this.pa, sum) == 0){
         this.pa <- this.sp
      }
      pa.stack[[i]] <- this.pa

      # Dropping out layers that have already been assigned
      temp.cluster <- dropLayer(temp.cluster, which.max(overlaps))
   }
   pa.stack <- dropLayer(pa.stack, 1)

   # Convert to 1/0
   pa.stack <- pa.stack > 0

   if(plot == TRUE){
      raster::plot(pa.stack)
   }

   print("sympatrifying...")

   # Here's where we're generating sympatry.  In the original sim I had this set up
   # to produce a 1 everywhere that was present in the PA raster, the true suitability
   # in any grid cell that was within the buffer but NOT present in the PA raster,
   # and zero everywhere else.  I can't recall why I did it that way though. I think
   # for simulation purposes it might be better to have the suitability at every grid cell
   # in the range and 0 elsewhere

   for(i in names(breadth)){
      temp.pa <- pa.stack[[i]]
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

   return(output)
}
