#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double fa2_distance(NumericVector& n1, NumericVector& n2)
{
  double dis = 0.0;
  for(int i=0;i<n1.size();i++)
    dis += pow(n1[i] - n2[i],2.0);
    
  return sqrt(dis);
}

/* Get ids */
// [[Rcpp::export]]
IntegerVector fa2_getids(IntegerMatrix& G)
{
  int N = G.nrow();
  int M, Gi;
  std::vector<int> ids;
  ids.push_back(G(0,0));
  for(int i=0;i<N;i++)
    for(int j=0;j<2;j++)
    {
      M  = ids.size();
      Gi = G(i,j);
      for(int k=0;k<M;k++)
      {
        if (Gi==ids[k]) break;
        if ((k+1)==M) ids.push_back(Gi);
      }
    }
    
  return wrap(ids);
}


/* Measures degree */
int fa2_degreei(IntegerMatrix& G, int id) 
{
  int d = 0;
  int N = G.nrow();
  int n1, n2;
  for(int i=0;i<N;i++)
  {
    n1 = G(i,0), n2 = G(i,1);
    if (n1 == n2) continue;
    d += ((n1 == id) | (n2 == id));
  }
    
  return d;  
}

// [[Rcpp::export]]
IntegerMatrix fa2_degree(IntegerMatrix& G) 
{
  
  IntegerVector ids = fa2_getids(G);
  IntegerMatrix d(ids.size(),2);
  d(_,0) = ids;
  
  for(int i=0;i<ids.size();i++)
    d(i,1) = fa2_degreei(G,ids[i]);
    
  return d;  
}

// [[Rcpp::export]]
IntegerVector fa2_degree_index(IntegerMatrix& G) 
{
  
  IntegerVector ids = fa2_getids(G);
  IntegerVector d(ids.size());
  
  for(int i=0;i<ids.size();i++)
    d[i] = fa2_degreei(G,ids[i]);
    
  return d;  
}

// [[Rcpp::export]]
NumericVector fa2_get_center(NumericMatrix& pos)
{
  int N = pos.nrow();
  int k = pos.ncol();
  NumericVector center(k);
  
  
  double dimrange[k][2];

  /* Initial minimum */
  for(int j=0;j<k;j++)
  {
    dimrange[j][0] = pos(0,j);
    dimrange[j][1] = pos(0,j);
  }
  
  /* Getting the range of the data */
  double tmp = 0.0;
  for(int i=1;i<N;i++)
    for(int j=0;j<k;j++)
    {
      tmp = pos(i,j);
      if (dimrange[j][0] > tmp) dimrange[j][0] = tmp;
      if (dimrange[j][1] < tmp) dimrange[j][1] = tmp;
    }
      
  /* Creating the center */
  for(int j=0;j<k;j++)
    center[j] = (dimrange[j][1] + dimrange[j][0])/2.0;
    
  return center;
}

/*
ATTRACTION REPULSION FORCES
*/

//' @param n1 Position of node 1
//' @param n2 Position of node 2
//' @param d1 Degree of node 1
//' @param d2 Degree of node 2
//' @param s1 Size of node 1
//' @param s2 Size of node 2
//' @param kr Repulsion constant 
//' @param kr2 
//' @param ka Attraction constant
double fa2_forcei(
  NumericVector n1, NumericVector n2,
  int d1, int d2,
  double s1=1.0, double s2=1.0,
  double kr=1.0, double kr2=100.0, double ka=1.0,
  bool nooverlap = false
  )
{
  double F = 0.0;
  
  if (nooverlap)
  {
    double dis = fa2_distance(n1,n2) - s1 - s2;
    
    if (dis > 0.0) F = kr*(double)((d1 + 1)*(d2 + 1))/dis + dis;
    else if (dis < 0.0) F = kr2*(double)((d1 + 1)*(d2 + 1));
  }
  else 
  {
    double dis = fa2_distance(n1,n2);
    F = kr*(double)((d1 + 1)*(d2 + 1))/dis + dis;
  }
  
  return(F);
}

//' Estimates the force of a system
// [[Rcpp::export]]
NumericMatrix fa2_force(
  NumericMatrix& pos,
  IntegerMatrix& G,
  NumericVector& size
  )
{
  int N    = G.rows();
  NumericMatrix F(N,N);
  
  // Computing degrees
  IntegerVector d = fa2_degree_index(G);
  double tmp=0.0;
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      if (i<j) 
      {
        tmp = fa2_forcei(pos(i,_),pos(j,_),d[i],d[j],size[i],size[j]);
        F(i,j) = tmp;
        F(j,i) = tmp;
      }
      
  return F;
}

/*** R
set.seed(1)

# Base graph
G    <- cbind(c(1,1:5),c(2,2,3,1,8,9))
pos  <- matrix(runif(nrow(G)*2),ncol=2)*10
size <- runif(nrow(G))
G
fa2_getids(G)

fa2_degree(G)

fa2_distance(pos[1,],pos[2,])

fa2_force(pos, G, size)
fa2_get_center(pos)
*/

