#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> dshl2(int nen, int nint){
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
            dshg[0][l] = 0.0;
            dshg[1][l] = 0.0;
        }

        if (nen == 3){
            dshg[0][l] = 1.0;
            dshg[1][l] = -2.0;
            dshg[2][l] = 1.0;
        }

        // TODO: Implementar casos quadráticos, cúbicos, etc. ...

    }

    return dshg;

}