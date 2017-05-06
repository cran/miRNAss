#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export(.calcThreshold)]]
double calcThreshold(const NumericVector& ye, const NumericVector& yt, int objective) {
    vector< pair<double,bool> > y;
    y.reserve(yt.size());

    double npos=0, nneg=0;
    unsigned int i;
    for(i=0; i<yt.size(); i++) {
        if(yt[i] != 0) {
            if(yt[i]>0) {
                y.push_back(pair<double,bool>(ye[i], true));
                npos++;
            }
            else {
                y.push_back(pair<double,bool>(ye[i], false));
                nneg++;
            }
        }
    }
    sort(y.begin(), y.end());

    double obj, best=0, TP=npos, TN=0;
    int lbest=-1, rbest=-1;
    if(objective == 0) { // Gm
        for(i=0; i<y.size()-1; i++) {
            if(y[i].second) TP--;
            else TN++;
            obj = TP / npos * TN / nneg;
            if(obj > best) {
                lbest=i; rbest=i+1;
                best = obj;
            }
            else if (obj == best)
                rbest=i;
        }
    } else if(objective==1) { // G
        for(i=0; i<y.size()-1; i++) {
            if(y[i].second) TP--;
            else TN++;
            obj = TP / npos * TP / (y.size()-i-1);
            if(obj > best) {
                lbest=i; rbest=i+1;
                best = obj;
            }
            else if(obj == best)
                rbest=i;
        }
    } else if(objective==2) { // F1
        for(i=0; i<y.size()-1; i++) {
            if(y[i].second) TP--;
            else TN++;
//             obj = 2 * TP / (2 * TP + (nneg - TN) + (npos - TP) );
            obj = 2 * TP / (TP - TN + y.size());
            if(obj > best) {
                lbest=i; rbest=i+1;
                best = obj;
            }
            else if(obj == best)
                rbest=i;
        }
    }
    return (y[lbest].first + y[rbest].first) / 2.0;
}

