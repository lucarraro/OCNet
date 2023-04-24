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

 #include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List continue_FD_SC(IntegerVector IndexHeadpixel, IntegerVector FD_to_SC, List SC_to_FD, IntegerVector downNode){
for (int ind{0}; ind<IndexHeadpixel.length(); ++ind){
int p = IndexHeadpixel[ind];
int pNew = p;
int k = 0;
IntegerVector sub_p;
while(k==0){
k = FD_to_SC[pNew-1];
if (k==0){
sub_p.push_back(pNew);
pNew = downNode[pNew-1];
}
}
IntegerVector foo (sub_p.length(), k);
FD_to_SC[sub_p-1] = foo;
IntegerVector tmp = SC_to_FD[k-1];
IntegerVector tmp2 (tmp.length() + sub_p.length());
int i=0;
for( ; i<tmp.size(); i++) tmp2[i] = tmp[i] ;
for( int j{0}; j<sub_p.size(); i++, j++) tmp2[i] = sub_p[j] ;

SC_to_FD[k-1] = tmp2;
}
List Lout = List::create(_["FD_to_SC"]=FD_to_SC, _["SC_to_FD"]=SC_to_FD);
return(Lout);
}