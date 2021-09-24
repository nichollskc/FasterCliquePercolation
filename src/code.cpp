#include <Rcpp.h>
using namespace Rcpp;

/*
 * Calculate clique communities roughly using the method outlined at the end
 * of Section 3 Farkas et al. The key idea is that to identify all the k-cliques
 * that neighbour a given clique, we can find all (k-1)-cliques which are contained
 * within the clique (subcliques) and then any k-clique containing one of these
 * subcliques will be a neighbour (if it exists in the network).
 */
class CliquePercolator {
  // The maximum node value we will consider i.e. when adding an extra node to
  // a subclique, we will try node=0,1,2,...,max_node_index
  int max_node_index;

  /* Set of unvisited cliques (i.e. they are present in the network but not
   * yet visited in the algorithm) */
  std::set<std::vector<int>> unvisited_cliques;

  /*
   * Keeps track of communities. Each entry in the vector corresponds to
   * a community, and the entry is the set of nodes that belong to that community.
   */
  std::vector<std::set<int>> community_membership;

  /*
   * Main function - returns a list of all the communities
   * Iteratively picks a new unvisited clique and then explores the community
   * containing this node (i.e. the connected component in the clique graph, which
   * has a node for each k-clique and nodes are connected if the k-cliques share
   * a (k-1)-clique.
   */
  public: std::vector<std::set<int>> calculate_community_membership() {
    // Will be the index into the community_membership vector
    int community_id = 0;

    // Pick a clique with which to start exploration
    std::set<std::vector<int>>::iterator next_clique = unvisited_cliques.begin();

    while (next_clique != unvisited_cliques.end()) {
      // Start a new community that will be spawned with this clique
      // The set of members starts off as an empty set
      std::set<int> empty;
      community_membership.push_back(empty);

      // Explore the community
      explore_community(*next_clique, community_id);

      // Pick the next clique and advance the community index
      community_id++;
      next_clique = unvisited_cliques.begin();
    }

    return(community_membership);
  }

  /*
   * Find all cliques which are connected to this clique (in the clique graph described
   * above) and mark all nodes contained in these cliques as a member of the
   * community with index community_id
   *
   * Mark all cliques in this community as visited by removing them from unvisited_cliques.
   */
  void explore_community(std::vector<int> clique, int community_id) {
    // To avoid recursion, maintain a stack (LIFO) of cliques whose neighbouring
    // cliques we will add to the community
    std::vector<std::vector<int>> cliques_to_check = {clique};

    // First mark this seed clique as visited
    unvisited_cliques.erase(unvisited_cliques.find(clique));

    while (!cliques_to_check.empty()) {
      // Remove a clique from the stack
      std::vector<int> current_clique = cliques_to_check.back();
      cliques_to_check.pop_back();

      // Mark all the nodes in this clique as belonging to the community we are exploring
      for (int j = 0; j < current_clique.size(); j++) {
        int node = current_clique[j];
        community_membership[community_id].insert(node);
      }

      // Find all adjacent cliques
      std::vector<std::vector<int>> adjacent = find_adjacent_unvisited_cliques(current_clique);
      for (int i = 0; i < adjacent.size(); i++) {
        // Add this adjacent clique to list of cliques to explore next
        cliques_to_check.push_back(adjacent[i]);
      }
    }
  }

  /*
   * Find all k-cliques which are adjacent to the clique in the clique graph. Mark
   * all such cliques as visited by removing them from unvisited_cliques and
   * return the list of adjacent cliques.
   */
  std::vector<std::vector<int>> find_adjacent_unvisited_cliques(std::vector<int> clique) {
    std::vector<std::vector<int>> adjacent_cliques;
    std::set<std::vector<int>>::iterator it;

    // Consider each (k-1) subclique of this clique
    int k = clique.size();
    for (int i = 0; i < k; i++) {
      // Create a copy and remove the ith element to create a (k-1)-clique
      // which is a subclique of clique
      std::vector<int> subclique = clique;
      subclique.erase(subclique.begin() + i);

      // Consider all k-cliques containing this (k-1) clique
      // All such k-cliques are neighbours of the original clique
      for (int j = 0; j <= max_node_index; j++) {
        // Create a copy of the subclique and add the potential new node
        std::vector<int> superclique = subclique;
        superclique.push_back(j);

        // Sort clique so that it is in the canonical representation
        sort(superclique.begin(), superclique.end());

        it = unvisited_cliques.find(superclique);
        if (it != unvisited_cliques.end()) {
          // This is a real, unvisited clique in the network
          adjacent_cliques.push_back(superclique);
          // Mark as visited by removing it
          unvisited_cliques.erase(it);
        }
      }
    }

    return(adjacent_cliques);
  }

  /*
   * Constructor - initialises necessary variables
   */
  CliquePercolator(std::vector<std::vector<int>> cliques, int max_node) {
    max_node_index = max_node;

    // Initialise empty list of communities
    std::vector<std::set<int>> no_communities;
    community_membership = no_communities;

    // Initialise unvisited cliques set to all cliques from list
    for (int i = 0; i < cliques.size(); i++) {
      unvisited_cliques.insert(cliques[i]);
    }
  }
};

// [[Rcpp::export]]
std::vector<std::set<int>> calculate_community_membership(std::vector<std::vector<int>> cliques, int max_node) {
  // Essentially just a wrapper - creates CliquePercolator object and then
  // calls its main function to find the communities.
  CliquePercolator cp = CliquePercolator(cliques, max_node);
  std::vector<std::set<int>> memberships = cp.calculate_community_membership();
  return(memberships);
}

// count matches between x and y, exploit fact they are pre-sorted
int count_matches(IntegerVector x,IntegerVector y) {
  int n=x.size();
  int i=0, j=0, matches=0;
  while(i < n && j < n) {
    if(x(i) < y(j))
      i++;
    else if(y(j) < x(i))
      j++;
    else if(x(i) == y(j)) {
      i++; j++; matches++;
    }
  }
  return(matches);
}

// [[Rcpp::export]]
IntegerMatrix create_edges_matrix_custom(List cliques) {
  int num_cliques=cliques.size();
  IntegerMatrix clique_adjacency(num_cliques, num_cliques);
  IntegerVector clique=cliques[0];
  int k=clique.size();
  Rcout << "n " << num_cliques << " k " << k << '\n';
  int num_edges=0;
  for(int i=0; i<num_cliques; i++) {
    IntegerVector clique_i=cliques[i];
    for(int j=(i+1); j<num_cliques; j++) {
      IntegerVector clique_j=cliques[j];
      int matches=count_matches(clique_i, clique_j);
      if(matches == (k-1)) {
        // Cliques i and j are "adjacent" so add this clique edge to the matrix
        clique_adjacency(i, j) = 1;
        clique_adjacency(j, i) = 1;
        num_edges++;
      }
    }
  }
  Rcout << " edges found " << num_edges << '\n';
  return( clique_adjacency );
}

// [[Rcpp::export]]
IntegerMatrix create_edges_matrix_intersect(List cliques) {
  int n=cliques.size();
  int nr=n*(n-1)/2;
  // if(n > 10000)
  //   nr= (int) sqrt( (double)nr ); // risky, but may save memory headaches
  IntegerMatrix edges(2,nr);
  IntegerVector clique=cliques[0];
  int k=clique.size();
  Rcout << "n " << n << " nr " << nr << " k " << k << '\n';
  int m=0;
  for(int i=0; i<n; i++) {
    IntegerVector clique_i=cliques[i];
    for(int j=(i+1); j<n; j++) {
      IntegerVector clique_j=cliques[j];
      IntegerVector ij_int = intersect( clique_i, clique_j );
      if(ij_int.size() == (k-1)) {
        // Rcout << m << ' ' << i << ' ' << j << '\n';
        edges(0,m)=i+1;
        edges(1,m)=j+1;
        m++;
      }
    }
  }
  Rcout << " edges found " << m << '\n';
  IntegerMatrix trimmed(2,0);
  if(m!=0)
  //   trimmed = edges( 0, Range(0,1) );
  // else
    trimmed = edges( Range(0,1), Range(0,(m-1)) );
  return( trimmed );
}

// [[Rcpp::export]]
IntegerMatrix create_edges_matrix_setdiff(List cliques) {
  int n=cliques.size();
  int nr=n*(n-1)/2;
  // if(n > 10000)
  //   nr= (int) sqrt( (double)nr ); // risky, but may save memory headaches
  IntegerMatrix edges(2,nr);
  IntegerVector clique=cliques[0];
  int k=clique.size();
  Rcout << "n " << n << " nr " << nr << " k " << k << '\n';
  int m=0;
  for(int i=0; i<n; i++) {
    IntegerVector clique_i=cliques[i];
    for(int j=(i+1); j<n; j++) {
      IntegerVector clique_j=cliques[j];
      IntegerVector ij_int = setdiff( clique_i, clique_j );
      if(ij_int.size() == 1) {
        // Rcout << m << ' ' << i << ' ' << j << '\n';
        edges(0,m)=i+1;
        edges(1,m)=j+1;
        m++;
      }
    }
  }
  Rcout << " edges found " << m << '\n';
  IntegerMatrix trimmed(2,0);
  if(m!=0)
  //   trimmed = edges( 0, Range(0,1) );
  // else
    trimmed = edges( Range(0,1), Range(0,(m-1)) );
  return( trimmed );
}
