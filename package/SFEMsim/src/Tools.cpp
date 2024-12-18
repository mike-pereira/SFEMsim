// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

using namespace Rcpp;

typedef Eigen::Triplet<double> T;

/*
 * Compute identity matrix
 */
Eigen::SparseMatrix<double> SparseId(int n){
  std::vector<T> tripletList;
  tripletList.reserve(n);
  for(int i=0; i<n;++i){
    tripletList.push_back(T(i,i,1));
  }
  Eigen::SparseMatrix<double> mat(n,n);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return mat;
}

/*
 * Compute the minimum value of a NumericVector
 */
double vecmin(NumericVector x) {
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  return *it;
}


// [[Rcpp::export]]
std::vector<int> argSort(const Eigen::ArrayXd& array) {
  // Create a vector of indices
  std::vector<int> indices(array.size());
  for (int i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }
  
  // Sort the indices based on the values in the Eigen array
  std::sort(indices.begin(), indices.end(),
            [&array](int a, int b) { return array[a] < array[b]; });
  
  return indices;
}



