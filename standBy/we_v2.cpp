#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<double> we(int nint){
    vector<double> we(nint);
    
    if(nint == 2){
        we[0] = 1.0;
        we[1] = 1.0;
    }

    if(nint == 3){
        we[0] = 0.8888888888888888;
        we[1] = 0.5555555555555556;
        we[2] = 0.5555555555555556;
    }

    if(nint == 4){
        we[0] = 0.6521451548625461;
        we[1] = 0.6521451548625461;
        we[2] = 0.3478548451374538;
        we[3] = 0.3478548451374538;
    }

    if(nint == 5){
        we[0] = 0.5688888888888889;
        we[1] = 0.4786286704993665;
        we[2] = 0.4786286704993665;
        we[3] = 0.2369268850561891;
        we[4] = 0.2369268850561891;
    }

    return we;

}