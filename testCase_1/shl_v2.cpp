#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> shl(int nen, int nint){
    vector<double> pt(nint);
    vector<double> w(nint);

    vector<vector<double>> shg(nint, vector<double>(nen));
    
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
            shg[0][l] = (1.0-t)/2.0;
            shg[1][l] = (1.0+t)/2.0;
        }

        if (nen == 3){
            shg[0][l] = 0.5*t*(t-1);
            shg[1][l] = -1.0*(t-1)*(t+1);
            shg[2][l] = 0.5*t*(t+1);
        }

        if (nen == 4){
            shg[0][l] = (-9.0/16.0)*(t+1.0/3.0)*(t-1.0/3.0)*(t-1.0);
            shg[1][l] = (27.0/16.0)*(t+1.0)*(t-1.0/3.0)*(t-1.0);
            shg[2][l] = (-27.0/16.0)*(t+1.0)*(t+1.0/3.0)*(t-1.0);
            shg[3][l] = (9.0/16.0)*(t+1.0)*(t+1.0/3.0)*(t-1.0/3.0);
        }

        if (nen == 5){
            shg[0][l] = (2.0/3.0)*(t+1.0/2.0)*t*(t-1.0/2.0)*(t-1.0);
            shg[1][l] = (-8.0/3.0)*(t+1.0)*t*(t-1.0/2.0)*(t-1.0);
            shg[2][l] = 4.0*(t+1.0)*(t+1.0/2.0)*(t-1.0/2.0)*(t-1.0);
            shg[3][l] = (-8.0/3.0)*(t+1.0)*t*(t+1.0/2.0)*(t-1.0);
            shg[4][l] = (2.0/3.0)*(t+1.0)*(t+1.0/2.0)*t*(t-1.0/2.0);
        }

        // TODO: Implementar casos quadráticos, cúbicos, etc. ...

    }

    return shg;

}