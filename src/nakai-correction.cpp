#include <Rcpp.h>

#ifndef PI 
#define PI 3.14159265
#endif

double sinerr(double x, double wd) {
// Sine correction function phi_sr(alpha, gamma)
// Eqs. (10), (13), and (14) of Nakai and Shimoyama (submitted)

    double a1, a2, a3, a4, a5;
    double b1, b2, b3, b4, b5;
    double f_aoa, x_orig, wd_orig;
    
    a1 = -3.19818998552857E-10;
    a2 = -2.69824417931343E-8;
    a3 = 4.16728613218081E-6;
    a4 = 4.85252964763967E-4;
    a5 = 1.67354200080193E-2;
    b1 = 5.92731123831391E-10;
    b2 = 1.44129103378194E-7;
    b3 = 1.20670183305798E-5;
    b4 = 3.92584527104954E-4;
    b5 = 3.82901759130896E-3;
    
    x_orig = x;
    wd_orig = wd; 
    
    if (x > 0) {
        x *= -1;
        wd += 180;
    }
    
    f_aoa = a1 * pow(x, 5) + a2 * pow(x, 4) + a3 * pow(x, 3) + a4 * pow(x, 2) + a5 * x + 1;
    f_aoa -= sin(3 * wd * PI/180) * (b1 * pow(x, 5) + b2 * pow(x, 4) + b3 * pow(x, 3) + b4 * pow(x, 2) + b5 * x);
    
    x = x_orig;
    wd = wd_orig;
    
    return f_aoa;
}

double coserr(double x, double wd) {
// Cosine correction function phi_cr(alpha, gamma)
// Eqs. (11), (12), (15), and (16) of Nakai and Shimoyama (submitted)

    double c1, c2, c3, c4, c5;
    double d1, d2, d3, d4, d5;
    double f_aoa, x_orig, wd_orig;
    
    c1 = -1.20804470033571E-9;
    c2 = -1.58051314507891E-7;
    c3 = -4.95504975706944E-6;
    c4 = 1.60799801968464E-5;
    c5 = 1.28143810766839E-3;
    d1 = 2.2715401644872E-9;
    d2 = 3.85646200219364E-7;
    d3 = 2.03402753902096E-5;
    d4 = 3.94248403622007E-4;
    d5 = 9.18428193641156E-4;

    x_orig = x;
    wd_orig = wd; 
    
    if (x > 0) {
        x *= -1;
        wd += 180;
    }
    
    if (x < -70)	x = -70;
    
    f_aoa = c1 * pow(x, 5) + c2 * pow(x, 4) + c3 * pow(x, 3) + c4 * pow(x, 2) + c5 * x + 1;
    f_aoa += sin(3 * wd * PI/180) * (d1 * pow(x, 5) + d2 * pow(x, 4) + d3 * pow(x, 3) + d4 * pow(x, 2) + d5 * x);
    
    x = x_orig;
    wd = wd_orig;
    
    return f_aoa;
}

double gx(double x, double wd, double a) {
// Equation to solve
    
    return atan(a * coserr(x, wd) / sinerr(x, wd)) * 180 / PI;
    
}

double Steffensen(double x, double wd, double a) {
// Appendix A of Nakai et al.(2006)

    double x0, x1, x2, x3;
    double key;
    
    x0 = x;
    while (1) {
        x1 = gx(x0, wd, a);
        x2 = gx(x1, wd, a);
        x3 = x2;	// output value if `break' out of the loop
        key = x2 - 2 * x1 + x0;
        if (fabs(key) < 0.01) break;
        x3 = x0 - pow(x1 - x0, 2) / key;
        x0 = x3;
    }
    
    return x3;
}

// [[Rcpp::export]]
Rcpp::List nakai_correction_2012(std::vector<double> &uIn, std::vector<double> &vIn, std::vector<double> &wIn) {
    const int N = uIn.size();
    Rcpp::NumericVector Uout(N);
    Rcpp::NumericVector Vout(N);
    Rcpp::NumericVector Wout(N);
    double aoa, sin_err, cos_err, wd, ws;
    int i;	
    for(int i = 0; i < N; i++){
        if (wIn[i] == 0){
            aoa = 0.0;
        } else {
            ws = sqrt(pow(uIn[i], 2) + pow(vIn[i], 2));
            
            if (ws == 0) {
                if (wIn[i] >= 0) aoa = 90;
                if (wIn[i] < 0)  aoa = -90;
            } else {
                aoa = atan(wIn[i]/ws) * 180 / PI;
            }
            
            // Wind direction
            if (ws == 0) {
                wd = 0;
            } else if (vIn[i] >= 0) {
                wd = 180 - acos(uIn[i] / ws) * 180 / PI;
            } else {
                wd = 180 + acos(uIn[i] / ws) * 180 / PI;
            }
            
            // Steffensen's method --- find out true AoA
            if (ws != 0) aoa = Steffensen(aoa, wd, wIn[i]/ws);
        }
        // sine error calculation
        sin_err = sinerr(aoa, wd);	
        // cosine error calculation
        cos_err = coserr(aoa, wd);
        Uout[i] = uIn[i] / cos_err;
        Vout[i] = vIn[i] / cos_err;
        Wout[i] = wIn[i] / sin_err;
    }
    return Rcpp::List::create(
        Rcpp::_["uCorr"] = Uout,
        Rcpp::_["vCorr"] = Vout,
        Rcpp::_["wCorr"] = Wout);
}

