
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List NN_OCN(int dimX, int dimY, bool periodicBoundaries, IntegerMatrix movement){
  List NeighbouringNodes(dimX*dimY);
  IntegerVector mov_r = movement(0, _);
  IntegerVector mov_c = movement(1, _);
  int k {0};
  IntegerVector v0 (8, 0);
  IntegerVector vX (8, dimX);
  IntegerVector vY (8, dimY);
  for (int cc{0}; cc<dimX; ++cc){
    for (int rr{0}; rr<dimY; ++rr){
      IntegerVector neigh_r(8, rr+1);
      IntegerVector neigh_c(8, cc+1);
      neigh_r = neigh_r + mov_r;
      neigh_c = neigh_c + mov_c;
      if (periodicBoundaries){
        for (int idr{0}; idr<neigh_r.length(); ++idr){
          if (neigh_r[idr]==0) neigh_r[idr]=dimX;
          if (neigh_r[idr]>dimY) neigh_r[idr]=1;
        }
        for (int idc{0}; idc<neigh_c.length(); ++idc){
          if (neigh_c[idc]==0) neigh_c[idc]=dimY;
          if (neigh_c[idc]>dimX) neigh_c[idc]=1;
        }
      }
      LogicalVector NAB = (neigh_r > v0) & (neigh_r <= vY) & (neigh_c > v0) & (neigh_c <= vX);
      neigh_c = (neigh_c - 1)*dimY;
      NeighbouringNodes[k] =  neigh_r[NAB] + neigh_c[NAB];
      k++;
    }
  }
  return(NeighbouringNodes);
}



#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List NN_river(int dimX, int dimY, bool periodicBoundaries, IntegerMatrix movement, IntegerVector toDEM, int nNodes){
  List NeighbouringNodes(dimX*dimY);
  IntegerVector mov_r = movement(0, _);
  IntegerVector mov_c = movement(1, _);
  IntegerVector v0 (8, 0);
  IntegerVector vX (8, dimX);
  IntegerVector vY (8, dimY);
  for (int ind{0}; ind<nNodes; ++ind){
    int nodeDEM = toDEM[ind];
    int cc = nodeDEM % dimX;
    if (cc==0) cc = dimX;
    int rr = (nodeDEM - cc)/dimX + 1;
    IntegerVector neigh_r(8, rr);
    IntegerVector neigh_c(8, cc);
    neigh_r = neigh_r + mov_r;
    neigh_c = neigh_c + mov_c;
    if (periodicBoundaries){
      for (int idr{0}; idr<neigh_r.length(); ++idr){
        if (neigh_r[idr]==0) neigh_r[idr]=dimX;
        if (neigh_r[idr]>dimY) neigh_r[idr]=1;
      }
      for (int idc{0}; idc<neigh_c.length(); ++idc){
        if (neigh_c[idc]==0) neigh_c[idc]=dimY;
        if (neigh_c[idc]>dimX) neigh_c[idc]=1;
      }
    }
    LogicalVector NAB = (neigh_r > v0) & (neigh_r <= vY) & (neigh_c > v0) & (neigh_c <= vX);
    neigh_r = (neigh_r - 1)*dimX;
    NeighbouringNodes[nodeDEM-1] =  neigh_r[NAB] + neigh_c[NAB];
  }
  return(NeighbouringNodes);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List NN_FD(int nNodes, int dimX, int dimY, List NeighbouringNodes, IntegerVector toDEM){ //  
  List NeighbouringNodes_FD(nNodes);
  IntegerVector DEM_to_FD(dimX*dimY);
  IntegerVector foo = seq(1,nNodes);
  DEM_to_FD[toDEM-1] = foo;
  for (int ind{0};ind<nNodes;++ind){
    int indDEM = toDEM[ind];
    IntegerVector tmp = NeighbouringNodes[indDEM-1];
    tmp = tmp - 1;
    tmp = DEM_to_FD[tmp];
    LogicalVector foo = (tmp > 0);
    NeighbouringNodes_FD[ind] = tmp[foo];
  }
  return(NeighbouringNodes_FD);
}

