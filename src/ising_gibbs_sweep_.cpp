#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix ising_gibbs_sweep_(IntegerMatrix state, NumericVector proba_beta){
    RNGScope scope;
    int size = state.rows();
    int s;
    int itop, ibottom, jright, jleft;
    GetRNGstate();
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            s = 0;
            itop = (i+1) % size;
            ibottom = ((i + size - 1) % size);
            jright = (j+1) % size;
            jleft = (j + size - 1) % size;
            s += state(itop, j)
            + state(ibottom, j)
            + state(i, jright)
            + state(i, jleft);
            state(i,j) = 2*((runif(1))(0) < proba_beta((s+4)/2)) - 1;
        }
    }
    PutRNGstate();
    return state;
}
