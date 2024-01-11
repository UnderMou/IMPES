#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> dshl(int nen, int nint){
    vector<double> pt(nint);
    vector<double> w(nint);

    vector<vector<double>> dshg(nint, vector<double>(nen));
    
    if(nint == 2){
        pt[0] = -0.5773502691896257;
        pt[1] = 0.5773502691896257;
    }

    if(nint == 3){
        pt[0] = 0.0;
        pt[1] = -0.7745966692414834;
        pt[2] = 0.7745966692414834;
    }

    if(nint == 4){
        pt[0] = -0.3399810435848563;
        pt[1] = 0.3399810435848563;
        pt[2] = -0.8611363115940526;
        pt[3] = 0.8611363115940526;
    }

    if(nint == 5){
        pt[0] = 0.0;
        pt[1] = -0.5384693101056831;
        pt[2] = 0.5384693101056831;
        pt[3] = -0.9061798459386640;
        pt[4] = 0.9061798459386640;
    }

    double t;
    for (int l = 0; l < nint; l++){
        t = pt[l];

        if (nen == 2){
            dshg[0][l] = -1.0/2.0;
            dshg[1][l] = 1.0/2.0;
        }

        if (nen == 3){
            dshg[0][l] = t - 0.5;
            dshg[1][l] = -2.0*t;
            dshg[2][l] = t + 0.5;
        }

        if (nen == 4){
            dshg[0][l] = (-27.0*pow(t,2) + 18.0*t + 1)/16.0;        // 0.0625 + 1.125*t - 1.68758*pow(t,2);
            dshg[1][l] = (9.0*pow(t,2) - 2.0*t - 3.0)*(9.0/16.0);   // -1.6875 - 1.125*t + 5.0625*pow(t,2);
            dshg[2][l] = (9.0*pow(t,2) + 2.0*t - 3.0)*(-9.0/16.0);  // 1.6875 - 1.125*t - 5.0625*pow(t,2);
            dshg[3][l] = (27.0*pow(t,2) + 18.0*t - 1.0)/16.0;       // -0.0625 + 1.125*t + 1.68758*pow(t,2);
        }

        if (nen == 5){
            dshg[0][l] = (16.0*pow(t,3) - 12.0*pow(t,2) -2.0*t + 1.0)/6.0;          // 0.166667 - 0.333333*t - 2.0*pow(t,2) + 2.66667*pow(t,3);
            dshg[1][l] = (8.0*pow(t,3) - 3*pow(t,2) - 4.0*t + 1.0)*(-4.0/3.0);      // -1.33333 + 5.33333*t + 4.0*pow(t,2) - 10.6667*pow(t,3);
            dshg[2][l] = (2.0*t)*(8.0*pow(t,2) - 5.0);                              // -10.0*t + 16.0*pow(t,3);
            dshg[3][l] = (8.0*pow(t,3) + 3.0*pow(t,2) - 4.0*t - 1.0)*(-4.0/3.0);    // 1.33333 + 5.33333*t - 4.0*pow(t,2) - 10.6667*pow(t,3);
            dshg[4][l] = (16.0*pow(t,3) + 12.0*pow(t,2) - 2.0*t - 1.0)/6.0;         // -0.166667 - 0.333333*t + 2.0*pow(t,2) + 2.66667*pow(t,3);
        }

        // TODO: Implementar casos quadráticos, cúbicos, etc. ...

    }

    return dshg;

}