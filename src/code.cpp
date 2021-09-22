#include <Rcpp.h>
using namespace Rcpp;

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
