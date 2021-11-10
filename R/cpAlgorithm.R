#' Clique Percolation Community Detection
#'
#' Function for clique percolation community detection algorithms for weighted
#' and unweighted networks.
#'
#' @param W A qgraph object; see also \link[qgraph]{qgraph}
#' @param k Clique size (number of nodes that should form a clique)
#' @param method A string indicating the method to use 
#'   (\code{"unweighted"}, \code{"weighted"}, or \code{"weighted.CFinder"}); see Details 
#' @param I Intensity threshold for weighted networks
#' 
#' @return
#'   A list object with the following elements:
#'   \describe{
#'   \item{list.of.communities.numbers}{list of communities with numbers as identifiers 
#'         of nodes}
#'   \item{list.of.communities.labels}{list of communities with labels from qgraph object 
#'         as identifiers of nodes}
#'   \item{shared.nodes.numbers}{vector with all nodes that belong to multiple communities 
#'         with numbers as identifiers of nodes}
#'   \item{shared.nodes.labels}{vector with all nodes that belong to multiple communities 
#'         with labels from qgraph object as identifiers of nodes}
#'   \item{isolated.nodes.numbers}{vector with all nodes that belong to no community 
#'         with numbers as identifiers of nodes}
#'   \item{isolated.nodes.labels}{vector with all nodes that belong to no community 
#'         with labels from qgraph object as identifiers of nodes}
#'   \item{k}{user-specified \code{k}}
#'   \item{method}{user-specified method}
#'   \item{I}{user-specified \code{I} (if method was \code{"weighted"} 
#'         or \code{"weighted.CFinder"})}
#' }
#'
#' @details
#'   \code{method = "unweighted"} conducts clique percolation for unweighted networks as
#'   described in Palla et al. (2005). \code{method = "weighted"} conducts clique percolation
#'   for weighted graphs with inclusion of cliques if their Intensity is higher than the
#'   specified Intensity (\code{I}), which is the method described in Farkas et al. (2007).
#'   \code{method = "weighted.CFinder"} conducts clique percolation as in the CFinder program.
#'   The Intensity (\code{I}) threshold is applied twice, namely first to the Intensity of the
#'   cliques (as before) and then also to their \code{k-1} overlap with other cliques
#'   (e.g., in the case of \code{k = 3}, it is applied to the edge that two cliques share).
#' 
#'   For weighted networks, the absolute value of the edge weights is taken.
#'   Therefore, negative edges are treated like positive edges just like in the CFinder program.
#'   Thus, the Intensity threshold \code{I} can only be positive.
#'
#'   cpAlgorithm produces a solution for all networks, even if there are no communities 
#'   or communities have no overlap. The respective output is empty in such cases.
#'
#' @examples
#' ## Example for unweighted networks
#' 
#' # create qgraph object
#' W <- matrix(c(0,1,1,1,0,0,0,0,
#'               0,0,1,1,0,0,0,0,
#'               0,0,0,0,0,0,0,0,
#'               0,0,0,0,1,1,1,0,
#'               0,0,0,0,0,1,1,0,
#'               0,0,0,0,0,0,1,0,
#'               0,0,0,0,0,0,0,1,
#'               0,0,0,0,0,0,0,0), nrow = 8, ncol = 8, byrow = TRUE)
#' W <- Matrix::forceSymmetric(W)
#' W <- qgraph::qgraph(W)
#'
#' # run clique percolation for unweighted networks
#' results <- cpAlgorithm(W = W, k = 3, method = "unweighted")
#' 
#' ## Example for weighted networks
#' 
#' # create qgraph object
#' W <- matrix(c(0,1,1,1,0,0,0,0,
#'               0,0,1,1,0,0,0,0,
#'               0,0,0,0,0,0,0,0,
#'               0,0,0,0,1,1,1,0,
#'               0,0,0,0,0,1,1,0,
#'               0,0,0,0,0,0,1,0,
#'               0,0,0,0,0,0,0,1,
#'               0,0,0,0,0,0,0,0), nrow = 8, ncol = 8, byrow = TRUE)
#' set.seed(4186)
#' rand_w <- stats::rnorm(length(which(W == 1)), mean = 0.3, sd = 0.1)
#' W[which(W == 1)] <- rand_w
#' W <- Matrix::forceSymmetric(W)
#' W <- qgraph::qgraph(W)
#' 
#' # run clique percolation for weighted networks
#' results <- cpAlgorithm(W = W, k = 3, method = "weighted", I = 0.1)
#' 
#' @references
#' Farkas, I., Abel, D., Palla, G., & Vicsek, T. (2007). Weighted network modules.
#' \emph{New Journal of Physics, 9}, 180-180. http://doi.org/10.1088/1367-2630/9/6/180
#' 
#' Palla, G., Derenyi, I., Farkas, I., & Vicsek, T. (2005). Uncovering the overlapping community 
#' structure of complex networks in nature and society. \emph{Nature, 435}, 
#' 814-818. http://doi.org/10.1038/nature03607
#' 
#' @author Jens Lange, \email{lange.jens@@outlook.com}
#' 
#' @importFrom magrittr "%>%"
#' 
#' @export cpAlgorithm

cpAlgorithm <- function(W, k, method = c("unweighted","weighted","weighted.CFinder"), I){
  ###error message if W is not a qgraph object
  if (methods::is(W, "qgraph") == FALSE) {
    stop("W (network object) must be a qgraph object.")
  }
  
  #extract weights matrix
  Wmat <- qgraph::getWmat(W)
  labels <- as.vector(W$graphAttributes$Nodes$labels)

  return(cpAlgorithmRaw(Wmat, k, method, I, labels, all_k_cliques = NULL))
}

cpAlgorithmRaw <- function(Wmat, k, method = c("unweighted","weighted","weighted.CFinder"), I, labels = NULL, all_k_cliques = NULL) {
  ###error message if k is not larger than 2
  if (k < 3) {
    stop("k must be larger than 2, because this is the first reasonable clique size.")
  }
  ###error message if method is not "unweighted", "weighted", or "weighted.CFinder"
  if (method != "unweighted" & method != "weighted" & method != "weighted.CFinder") {
    stop("method must be 'unweighted', 'weighted', or 'weighted.CFinder' (depending on the network).")
  }
  
  ## ###function to determine list of cliques for unweighted networks
  ## unweighted <- function(W_unweighted){
    
  ##   #transform to igraph object to use specific functions
  ##   W_i_unweighted <- igraph::graph_from_adjacency_matrix(W_unweighted,
  ##                                                         mode = "undirected", weighted = NULL)
    
  ##   #extract cliques
  ##   cliques_unweighted <- igraph::cliques(W_i_unweighted, min = k, max = k) %>% lapply(as.vector)
    
  ##   #return cliques
  ##   return(cliques_unweighted)
  ## }
  
  
  #function to derive communities, shared nodes, and isolated nodes
  results_cp <- function(Wmat, cliques, labels){
    
    #loop to compare all cliques with each other
    #if cliques share k-1 nodes, a vector is created stating their indices
                                        #if there are not at least two cliques, there can be no edge; thus, the edge list is empty
    #CFinder applies the Intensity threshold twice, once for the cliques and once for the overlap of cliques
    #this is not stated in the paper, but to increase comparability, this is also implemented here for method = weighted.CFinder
    #each k-1 set of nodes coded by an edge is taken and this subnetwork is again tested against I
    #each subnetwork coded by an edge is excluded, when it does not exceed I
    #this is done only when there are edges
    ## if (method == "weighted.CFinder" & length(edges) > 0) {
    ##   intensity_weighted_edge_net <- list()
    ##   for (i in 1:length(edges)) {
    ##     edge_net_nodes <- Reduce(intersect, list(cliques[[edges[[i]][1]]],cliques[[edges[[i]][2]]]))
    ##     edge_net_W <- abs(W)[edge_net_nodes,edge_net_nodes]
    ##     weights_edge_net <- c(edge_net_W[upper.tri(edge_net_W)])
    ##     exponent_edge_net <- 2/(length(edge_net_nodes) * (length(edge_net_nodes) - 1))
    ##     intensity_weighted_edge_net[[i]] <- prod(weights_edge_net)^exponent_edge_net
    ##   }
    ##   edges_include_weighted <- lapply(intensity_weighted_edge_net, function(x) as.numeric(as.character(x)) > I)
    ##   edges <- edges[which(edges_include_weighted == TRUE)]
    ## }
    
    if (length(cliques) > 1) {
      communities = calculate_community_membership(cliques, nrow(Wmat))
    } else {
      #communities list is empty
      communities <- list() 
    }
    
    #create vector of nodes that do not belong to a community
    #if there are no isolated nodes, create empty variable
    isolated <- subset(1:nrow(Wmat),
                       subset = !(1:nrow(Wmat)%in%unique(unlist(communities))))
    if (length(isolated) == 0) {isolated <-  c()}
    
    #if there is more than one community...
    #get shared nodes
    #if there are no shared nodes, create empty variable
    if (length(communities) > 1) {
      shared <- list()
      count <- 1
      for (i in 1:(length(communities) - 1)) {
        for (j in (i + 1):length(communities)) {
          shared[[count]] <- intersect(communities[[i]],communities[[j]])
          count <- count + 1
        }
      }
      shared <- sort(unique(unlist(shared)))
      if (length(shared) == 0) {shared <- c()}
    }
    
    #if there are communities...
    #create list of communities with labels from qgraph object instead of node numbers
    if (length(communities) > 0) {
      communities_labels <- vector("list",length(communities))
      for (i in 1:length(communities)) {
        communities_labels[[i]] <- labels[communities[[i]]]
      }
    }
    #if there are no communities...
    if (length(communities) == 0) {
      communities_labels <- list()
    }
    
    #create vector of isolated nodes with labels from qgraph object instead of node numbers
    #if there are no isolated nodes, create empty variable
    isolated_labels <- labels[isolated]
    if (length(isolated_labels) == 0) {isolated_labels <- c()}
    
    #if there is more than one community...
    #create vector with shared nodes with labels from qgraph object instead of node numbers
    if (length(communities) > 1) {
      shared_labels <- labels[shared]
    }
    #if there is less than two communities
    if (length(communities) < 2) {
      shared <- c()
      shared_labels <- c()
    }
    #if there are no shared nodes, create empty variable
    if (length(shared_labels) == 0) {shared_labels <- c()}
    
    #return community lists, shared nodes vectors, and isolated nodes vectors
    return(list(communities,communities_labels,
                shared,shared_labels,
                isolated,isolated_labels))
  }
  
  #run corresponding functions for respective method
  
  ## if (method == "unweighted") {
  ##   cliques_unweighted <- unweighted(Wmat)
  ##   results_unweighted <- results_cp(Wmat, cliques_unweighted, labels)
  ##   names(results_unweighted) <- c("list.of.communities.numbers","list.of.communities.labels",
  ##                                  "shared.nodes.numbers","shared.nodes.labels",
  ##                                  "isolated.nodes.numbers","isolated.nodes.labels")
  ##   return(c(results_unweighted,list(k = k, method = method)))
  ## }
  
  if (method == "weighted" | method == "weighted.CFinder") {
    ###error message if I is not larger than zero
    if (I <= 0) {
      stop("Intensity (I) must be larger than zero.\nThis is because all edges are considered positive.\nComparing intensities of cliques or their overlap to zero or negative value therefore makes no sense.")
    }

    # First find all cliques (if not provided as an argument), with their associated intensities
    if (is.null(all_k_cliques)) {
        all_k_cliques <- calculate_all_clique_intensities(Wmat, k)
    }

    # Restrict to cliques above this threshold
    cliques_above_thresh <- threshold_cliques(all_k_cliques, I)$cliques

    # Run algorithm (i.e. find communities)
    results_weighted <- results_cp(Wmat, cliques_above_thresh, labels)
    names(results_weighted) <- c("list.of.communities.numbers","list.of.communities.labels",
                                 "shared.nodes.numbers","shared.nodes.labels",
                                 "isolated.nodes.numbers","isolated.nodes.labels")
    return(c(results_weighted,list(k = k, method = method, I = I)))
  }
  
}

threshold_cliques <- function(cliques_result, I) {
  include_clique <- lapply(cliques_result$intensity_weighted, function(x) as.numeric(as.character(x)) > I)
  cliques_weighted_thr <- cliques_result$cliques[which(include_clique == TRUE)]
  intensity_weighted_thr <- cliques_result$intensity_weighted[which(include_clique == TRUE)]
  return(list("intensity_weighted" = intensity_weighted_thr,
              "cliques" = cliques_weighted_thr))
}

##function to find all cliques and calculate their intensities
calculate_all_clique_intensities_raw <- function(W_weighted, k) {
  #take absolute Value of weights matrix
  #deals with negative edges such that they are simply considered like positive edges
  W_weighted <- abs(W_weighted)

  #transform to igraph object to use specific functions
  W_i_weighted <- igraph::graph_from_adjacency_matrix(W_weighted,
                                                      mode = "undirected", weighted = TRUE)

  #extract cliques
  cliques_weighted <- igraph::cliques(W_i_weighted, min = k, max = k) %>% lapply(as.vector)

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
  }

  return(list("intensity_weighted" = intensity_weighted,
              "cliques" = cliques_weighted))

}

# Use memoization to cache the results of this function call. Thus if the function is called
# twice in a session with the same parameters, the second call will be almost instant.
calculate_all_clique_intensities <- R.cache::addMemoization(calculate_all_clique_intensities_raw)
