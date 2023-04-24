#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List WSC(int nNodes, List SC_to_FD, IntegerVector FD_to_SC, List NeighbouringNodes){
IntegerVector indices1;
IntegerVector indices2;
 for (int i{0}; i<nNodes; ++i)
  {
  IntegerVector set = SC_to_FD[i];
  IntegerVector nodes;

  for (int s{0}; s<set.length(); ++s)
  {
  IntegerVector nn = NeighbouringNodes[set[s]-1];
  nn = nn - 1;
  nn = FD_to_SC[nn]; 
  nodes = union_(nodes, nn);
  }
  IntegerVector tmp {i+1};
  IntegerVector neighSubcatch = setdiff(nodes, tmp);
  for (int kk{0}; kk<neighSubcatch.length(); ++kk){
  indices1.push_back(i+1);
  indices2.push_back(neighSubcatch[kk]);
  }
  }
  List Lout = List::create(_["ind1"] = indices1, _["ind2"] = indices2);
  return(Lout);
}