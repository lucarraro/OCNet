#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List paths_cpp(S4 OCN, IntegerVector whichNodes, String str = "RN", bool includePaths = false, bool includeDownstreamNode = false, bool includeUnconnectedPaths = false){
  
  //
  List L = OCN.slot(str);
  int nNodes = L["nNodes"];
  LogicalVector whichNodes_bool (nNodes);
  for (int i{0}; i < whichNodes.length(); ++i)
  {
	  int ind = whichNodes[i];
	  ind = ind - 1;
	  whichNodes_bool[ind] = true;
  }
  NumericVector downNode = L["downNode"];
  NumericVector foo = NumericVector::create(1000,0.1*nNodes);
  foo = ceiling(foo);
  int siz = nNodes*max(foo);
  IntegerVector set_col ( siz );
  IntegerVector set_row ( siz );
  NumericVector set_values ( siz );
  NumericMatrix downstreamLengthUnc;    // initialize in any case
  if (includeUnconnectedPaths)
  {
    NumericMatrix foo (nNodes); 
    downstreamLengthUnc = foo;
  }
  
  NumericVector leng = L["leng"];
  List upstream = L["upstream"];
  IntegerVector outlet = L["outlet"];
  int k = 0;
  
  List downstreamPath(nNodes); // initialize even when includePaths = false, otherwise the compiler complains
  if (includePaths) // initialize downstreamPath
  {
    for (int i{0}; i < nNodes; ++i)
    {
      List foo (nNodes);
      downstreamPath[i] = foo;
    }
  }
  
  for (int ind1{0}; ind1<whichNodes.length(); ++ind1)
  {
	  int i = whichNodes[ind1];
	  i = i - 1;
    if (includePaths)
      as<List>(downstreamPath[i])[i] = i+1;
    if (includeDownstreamNode)
    {
      set_row(k) = i;
      set_col(k) = i;
      set_values(k) = leng(i);
      k++;
    }
    IntegerVector path = IntegerVector::create(i+1);
    int j = i;
    int node_j = j+1;
	double sl = 0;
    while (!(std::find(outlet.begin(), outlet.end(), node_j)!=outlet.end())) //(!(node_j == outlet))
    {
      node_j = downNode(j);
      j = node_j-1;
      path.push_back(node_j);
	  NumericVector tmp = leng[path-1];
	  sl = sum(tmp);
	  if (whichNodes_bool[j])
	  {
		if (includePaths)
			as<List>(downstreamPath[i])[j] = path;
		set_row(k) = i;
		set_col(k) = j;
		if (includeDownstreamNode)
			set_values(k) = sl;
		else
			set_values(k) = sl - leng(j);
		k++;
      }
	  
      if (includeUnconnectedPaths)
      {
        IntegerVector tmp1 = upstream[j];
		tmp1 = intersect(tmp1, whichNodes);
        int ind = path[path.length()-2] -1;
        IntegerVector tmp2 = upstream[ind];
        tmp2.push_back(node_j);
        IntegerVector ups = setdiff(tmp1,tmp2);
        for (int ind_u{0}; ind_u<ups.length(); ++ind_u)
        {
          int node_u = ups[ind_u];
          int u = node_u - 1;
          if (includeDownstreamNode)
            downstreamLengthUnc(i, u) = sl;
          else
            downstreamLengthUnc(i, u) = sl - leng(j);
        }
      }
      
    }
  }
  
  set_values = set_values[ Range(0, k-1) ];
  set_row = set_row[ Range(0, k-1) ] + 1;
  set_col = set_col[ Range(0, k-1) ] + 1;
  List Lout = List::create( _["i"] = set_row, _["j"] = set_col, _["values"] = set_values);
  if (includePaths)
    Lout.push_back(downstreamPath, "downstreamPath");
  if (includeUnconnectedPaths)
    Lout.push_back(downstreamLengthUnc, "downstreamLengthUnc");
  
  return(Lout);
}
