###function to determine list of cliques for weighted networks
old_weighted <- function(W_weighted, I, k){

  #take absolute Value of weights matrix
  #deals with negative edges such that they are simply considered like positive edges
  W_weighted <- abs(W_weighted)

  #transform to igraph object to use specific functions
  W_i_weighted <- igraph::graph_from_adjacency_matrix(W_weighted,
                                                      mode = "undirected", weighted = TRUE)

  #extract cliques
  cliques_weighted <- igraph::cliques(W_i_weighted, min = k, max = k) %>% lapply(as.vector)

  #if there are cliques...(if there are no cliques, nothing needs to be done because cliques_weighted is empty)
  #check whether cliques exceed Intensity
  #first step is to create a list of intensity with each vector being the intensity of corresponding clique
  #this is achieved by extracting the weights for every pair of nodes in each clique from weights matrix
  #the vector with the weights is then used to calculate the intensity of the respective clique
  #second step is to create a list that has TRUE at corresponding position if clique intensity exceeds set Intensity
  #finally, select only cliques that should be included
  ## intensity_weighted <- list()
  intensity_weighted <- vector("list",length(cliques_weighted))
  if (length(cliques_weighted) > 0) {
    for (i in 1:length(cliques_weighted)) {
      ## weights <- c()
      weights <- numeric( (length(cliques_weighted[[i]])-1) * (length(cliques_weighted[[i]])) / 2 )
      m <- 1
      for (j in 1:(length(cliques_weighted[[i]]) - 1)) {
        for (k in (j+1):length(cliques_weighted[[i]])) {
          weights[m] <- W_weighted[cliques_weighted[[i]][j],cliques_weighted[[i]][k]]
          m <- m + 1
        }
      }
      exponent <- 2/(length(cliques_weighted[[i]]) * (length(cliques_weighted[[i]]) - 1))
      intensity_weighted[[i]] <- prod(weights)^exponent
    }
    cliques_include_weighted <- lapply(intensity_weighted, function(x) as.numeric(as.character(x)) > I)
    cliques_weighted <- cliques_weighted[which(cliques_include_weighted == TRUE)]
  }

  #return cliques
  return(list(cliques=cliques_weighted,
              all_intensities=intensity_weighted))

}

adj=structure(c(0, 0.01624581, 0.01541721, 0.01547624, 0.01624581,
                0, 0.79870296, 0.81206211, 0.01541721, 0.79870296, 0, 0.93482651,
                0.01547624, 0.81206211, 0.93482651, 0), .Dim = c(4L, 4L), .Dimnames = list(
                    c("A", "B", "C", "D"), c("A", "B", "C", "D")))
W <- qgraph::qgraph(adj, DoNotPlot=TRUE)
graphs = list(single_clique=W)

test_that("same cliques", {
  Wmat = qgraph::getWmat(graphs$single_clique)
  k = 3
  I = 0.5
  old_result = old_weighted(Wmat, I, k)
  new_cliques_with_intensities = calculate_all_clique_intensities(Wmat, k)
  thresholded_cliques = threshold_cliques(new_cliques_with_intensities, I)$cliques

  # At the very least should find same number of cliques
  expect_equal(length(old_result$all_intensities), length(new_cliques_with_intensities$intensity_weighted))
  # Since we haven't changed the code that much, actually expect same order of cliques so can check intensities simply
  expect_equal(old_result$all_intensities, new_cliques_with_intensities$intensity_weighted)

  # Check end lot of cliques after thresholding
  expect_equal(old_result$cliques, thresholded_cliques)
})
