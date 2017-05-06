#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export(.calcDistance)]]
NumericVector calcDistance(const IntegerVector& p,
                           const IntegerVector& idx,
                           const NumericVector& x,
                           const NumericVector& y) {
    NumericVector d(y.begin(), y.end());
    double val;
    int i, j;
    bool haveChange=true;

    while(haveChange) {
        haveChange=false;
        for(i=0; i<d.size(); i++) {
            for(j=p[i]; j<p[i+1]; j++) {
                val = d[idx[j]] * x[j];
                if(val > d[i]) {
                    haveChange=true;
                    d[i] = val;
                }
            }
        }
    }
    return d;
}
