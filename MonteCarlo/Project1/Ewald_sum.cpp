/**
 * @file Ewald_sum.cpp
 * @author Seongsu, Kim (you@domain.com)
 * @brief Ewald method header file that is made for caculation of Ising model
 * @version 0.1
 * @date 2022-03-22
 * 
 * @copyright Copyright (c) 2022
 * @see
 * Gamma : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/sf_gamma/tgamma.html
 * Incomplete Gamma : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
 * Exponent integral : https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/expint/expint_i.html
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <iostream>

using namespace std;

double pi = 3.1415926539;

// Periodic boundary condition
class ewald_ND{
    
    public:
        int d; // Dimension of the Ewald model.. but it looks useless
        double alpha; // alpha is "d + sigma" (overall inverse power)
        double kappa;
        
        vector<int> L; // Lattice size of each edge L[0] = Lx, L[1] = Ly

        ewald_ND(){};
        ewald_ND(int kd, vector<int> kL, double kalpha);

        int xof(int i){return i%L[0];}
        int yof(int i){return i/L[0];}

        // For general dimension Ndof = N dimensional of, Ndofh = N dimensional of helper
        int Ndof(int i, int d){return Ndofh(i,d)%L[d];}
        int Ndofh(int i,int d){
            if(d == 0) return i;
            else return Ndofh(i,d-1)/L[d-1];
        }

        // 2D function
        double dist_2pt(double x1, double y1, double x2, double y2){return pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2),0.5);};
        double dist_ij(int i, int j){return dist_2pt(xof(i),yof(i),xof(j),yof(j));};
        double dist_ij_min(int i, int j){
            double diffx = min(abs(xof(i)-xof(j)),L[0]-abs(xof(i)-xof(j)));
            // cout << diffx;
            double diffy = min(abs(yof(i)-yof(j)),L[1]-abs(yof(i)-yof(j)));
            // cout << diffy;
            return pow(diffx*diffx + diffy*diffy,0.5);
        };
        double vec_mul(double x1, double y1, double x2, double y2){return x1*x2 + y1* y2;};
        // int cache_point_2D(int i, int j){return (N)*i-i*(i+1)/2+(j-i-1);};
        int cache_point_2D(int i, int j){
            int nx = abs(xof(i)-xof(j)); nx = 2*nx > L[0] ? L[0] - nx : nx;
            int ny = abs(yof(i)-yof(j)); ny = 2*ny > L[0] ? L[0] - ny : ny;
            return nx+ny*((L[0])/2+1);
            // return (abs(abs(xof(i)-xof(j))-(L[0])/2))+(abs(abs(yof(i)-yof(j))-(L[1])/2))*(L[0]/2);
        };

        double pi_1(int i, int j, double kappa);
        double pi_2(int i, int j, double kappa);
        double pi_ij(int i, int j);

        double pi_1_1D(int i, int j, double kappa);
        double pi_2_1D(int i, int j, double kappa);
        double pi_ij_1D(int i, int j, bool consider_parallel = true);
        double periodic_sum_comparison_1D(int i, int j, int idx);

        double gam(double z);
        double gaminc_up(double z, double alpha);

        vector<double> dist_cache_2D; // O(N^4) 메모리 크기는 좀 힘들지 않을까-> v.max_size = 1e19 나오는디..?
                                   // N= 1000까지도 가능할듯? -> 실제로 해보려니까 뻗어버리네..
                                   // dist_cache_2D 크기는 n/4만 해도 가능 한 듯
        vector<double> dist_cache_1D;
};

ewald_ND::ewald_ND(int kd, vector<int> kL, double kalpha): d(kd){ //kL: x,y,z
    alpha = kalpha;
    L = kL;
    kappa = 2/(double)L[0]; // Multi dimension 에서는 다르게 해야할 듯..?

    // Memoization
    unsigned long long size = ((L[0])/2+1)*((L[1])/2+1);
    // cout << dist_cache_2D.max_size() << ' ' << N << ' ' << size << '\n';
    dist_cache_2D = vector<double>(size,-1);
    dist_cache_1D = vector<double>(L[1],-1);
}

// Gamma function
double ewald_ND::gam(double z){
    return boost::math::tgamma(z);
}

double ewald_ND::gaminc_up(double alpha, double z){ // not normalized
    if(alpha > 0) return boost::math::tgamma(alpha,z);
    else if(alpha == 0) return -(boost::math::expint(-z)); // if(alpha == 0) return -Ei(-z)
    else{ // negative alpha is not defined in cpp boost
        double sum_res = 0;
        double mul_res = 1;

        double s = alpha;
        int n = -floor(s);
        // cout << s << ' ' << n << ' ' << z << '\n';

        for(int i = 1; i <= n; i++){
            double temp = pow(z,s+i-1)*exp(-z);
            for(int j = 1; j <= i; j++) temp *= 1/(double)(s+j-1);
            sum_res += temp;

            mul_res *= 1/(s+i-1);
        }
        mul_res *= gaminc_up(s+n,z);
        // cout << -sum_res << ' ' << mul_res << '\n';
        return (-sum_res) + mul_res;
    }
}

double ewald_ND::pi_1(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 0; idx <= 5; idx++){
        for(int ii = 0; ii <= idx; ii++){ // nx, ny are nu_x and nu_y
        int nx = ii, ny = idx - ii; // if idx = 2, (nx, ny) = (0, 2) (1, 2) (2, 2)
            for(int kk = 0; kk < 2-(nx==0); kk++)
            for(int ll = 0; ll < 2-(ny==0); ll++){
                // dis = |r_i - r_j + L|
                double dis = dist_2pt(xof(i)-xof(j),yof(i)-yof(j),nx*L[0]*(1-2*kk),ny*L[1]*(1-2*ll));
                double inv_dist = pow(dis,-alpha);
                double temp = inv_dist*gaminc_up(alpha/2,(kappa*dis)*(kappa*dis));
                res += temp;
            }
        }        
    }
    return 1/gam(alpha/2)*res;
}

double ewald_ND::pi_2(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 1; idx <= 5; idx++){ // idx = 0 should be skipped? Maybe.. in thesis, h=0 should be vanished
        for(int ii = 0; ii <= idx; ii++){ // nx, ny are nu_x and nu_y
        double nx = (double)ii/L[0], ny = (double)(idx - ii)/L[1]; // if idx = 2, (nx, ny) = (0, 2) (1, 2) (2, 2)
            for(int kk = 0; kk < 2-(nx==0); kk++)
            for(int ll = 0; ll < 2-(ny==0); ll++){
                double hx = nx*(1-2*kk), hy = ny*(1-2*ll);
                double k_res = cos(2*pi*vec_mul(hx,hy,xof(i)-xof(j),yof(i)-yof(j)));
                double pi_k = pi*dist_2pt(hx,hy,0,0); 
                // k_res *= pow(0.5*pi_k,alpha-d);
                k_res *= 0.5*pow(pi_k,alpha-d);

                k_res *= gaminc_up(-0.5*(alpha-d),(pi_k/kappa)*(pi_k/kappa));
                // cout << idx << ' ' << hx << ' ' << hy << ' ' << k_res <<'\n';
                res += k_res;
            }
        }        
    }
    return 2*pow(pi,d/2)/gam(alpha/2)/(L[0]*L[1])*res;
}

double ewald_ND::pi_ij(int i,int j){
    if(i == j) return 0;
    // else if(i > j) swap(i,j); // i < j always true
    double &res = dist_cache_2D[cache_point_2D(i,j)];
    // cout << cache_point_2D(i,j) << '\n';
    if(res == -1){
        res = pi_1(i,j,kappa) + pi_2(i,j,kappa);
    }
    return res;
    // return pi_1(i,j,kappa) + pi_2(i,j,kappa);
}

// 만약에 parallel한 성분도 고려해야한다면 p1_1_1D, pi_2_2D 모두 다 pararell 방향의 mirror 이미지도 고려할 수 있게 수정이 필요함
double ewald_ND::pi_1_1D(int i, int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    for(int idx = 0; idx <= 5; idx++){
        int nx = idx;
        for(int kk = 0; kk < 2-(nx==0); kk++){ // This loop is used for convert sign of nx
            nx *= -1;
            double dis = dist_2pt(yof(i)-yof(j),0,-nx*L[1],0);
            // double dis = pow(abs(yof(j)-yof(i)+nx*L[1]),0.5);
            // cout << dis1 << " " << dis << '\n';
            double inv_dist = pow(dis,-alpha);
            double temp = inv_dist*gaminc_up(alpha/2,(kappa*dis)*(kappa*dis));
            res += temp; 
            // cout << idx <<" " << yof(i)-yof(j) <<  " " << nx << " " << temp << '\n';
        }
    }
    return 1/gam(alpha/2)*res;
}

double ewald_ND::pi_2_1D(int i,int j, double kappa){
    double res = 0;
    // max(idx) is 4 which is suggesetd at thesis. it makes cutoff 1e-15
    // But I set it as 5 for sure
    int d = 1;
    for(int idx = 1; idx <= 5; idx++){ // idx = 0 should be skipped? Maybe.. in thesis, h=0 should be vanished
        double nx = (double)idx/L[1];
        for(int kk = 0; kk < 2; kk++){ // This loop is used for convert sign of nx
            nx *= -1;
            double k_res = cos(2*pi*vec_mul(nx,0,yof(i)-yof(j),0));
            // double pi_k = pi*dist_2pt(nx,0,0,0);
            double pi_k = pi*abs(nx); 

            k_res *= 0.5*pow(pi_k,alpha-d);
            k_res *= gaminc_up(-0.5*(alpha-d),(pi_k/kappa)*(pi_k/kappa));
            // cout << idx << ' ' << hx << ' ' << hy << ' ' << k_res <<'\n';
            res += k_res;
            // cout << idx << " " << nx << " " << k_res << '\n';
        }
    }
    return 2*pow(pi,d/2)/gam(alpha/2)/L[1]*res;
}

double ewald_ND::pi_ij_1D(int i, int j, bool consider_parallel){
    if(i == j) return 0;
    if(consider_parallel && (xof(i) != xof(j))) return 0;
    double &res = dist_cache_1D[abs(yof(i)-yof(j))];
    if(res == -1){
        res = pi_1_1D(i,j,kappa) + pi_2_1D(i,j,kappa);

    }
    return res;
}

double ewald_ND::periodic_sum_comparison_1D(int i, int j, int idx){
    double res = 0;
    int dis = yof(j)-yof(i); // just difference
    res += pow(abs(dis),-alpha);
    for(int k = 1; k < idx; k++){
        res += pow(abs(dis+L[0]*k*1),-alpha);
        res += pow(abs(dis+L[0]*k*-1),-alpha);
    }
    return res;
}

// (1) Ewald sum number comparison generator
// int main(){
//     int lat = 16;
//     ewald_ND e2d(2,vector<int>({lat,lat}),3);
//     double kappa = 2/(double)lat;
//     // cout << e2d.kappa << '\n';
//     int cnt = 0;
//     for(int i = lat; i < lat*lat; i+=lat){
//         double t1 = e2d.pi_1_1D(0,i,kappa);
//         double t2 = e2d.pi_2_1D(0,i,kappa);
//         // cout << i << ' ' << e2d.cache_point_2D(0,i) << '\n'; 
//         cout << i << ' ' << e2d.periodic_sum_comparison_1D(0,i,30) << ' ' << t1 << " " << t2 << ' ' << e2d.pi_ij_1D(0,i) << ' ' << ((e2d.pi_ij_1D(0,i) - t1-t2) < 1e-20) << '\n';
//         // cout << i << ' ' << e2d.periodic_sum_comparison_1D(0,i,30) << ' ' << t1 << " " << t2 << ' ' << '\n';
//         // if(!((e2d.pi_ij_1D(0,i) - t1-t2)< 1e-19)) cnt++;
//     }
//     cout << cnt << '\n';
// }

// int main(){
//     ewald_ND e2d(2,vector<int>({256,256}),3);
//     double kappa = 2/(double)256;
//     cout << e2d.kappa << '\n';
//     int cnt = 0;
//     for(int i = 1; i < 256*256; i++){
//         double t1 = e2d.pi_1(0,i,kappa);
//         // // cout << t1 << '\n';
//         double t2 = e2d.pi_2(0,i,kappa);
//         // cout << t2 << '\n';
//         // cout << i << ' ' << t1 << ' ' << t2 << ' ' << t1 + t2 << '\n'; 
//         cout << i << ' ' << e2d.cache_point_2D(0,i) << '\n'; 
//         cout << i << ' ' << t1 + t2 << ' ' << e2d.pi_ij(0,i) << ' ' << ((e2d.pi_ij(0,i) - t1-t2) < 1e-20) << '\n';
//         if(!((e2d.pi_ij(0,i) - t1-t2)< 1e-19)) cnt++;
//         // cout << i << ' ' << e2d.pi_ij(0,i) << '\n'; 
//         // cout << abs(abs(e2d.xof(0)-e2d.xof(i))-e2d.L[0]/2) << '\n';
//         // cout << (e2d.L[0]/2-abs(abs(e2d.xof(0)-e2d.xof(i))-e2d.L[0]/2)) << ' ' << (e2d.L[1]/2-abs(abs(e2d.yof(0)-e2d.yof(i))-e2d.L[1]/2))*e2d.L[0] << '\n'; 
//     }
//     cout << cnt << '\n';
//     // cout << boost::math::tgamma(2,3) << '\n';
//     // cout << e2d.gaminc_up(2,3) << '\n';
// }

// 1) print_pi_ij
double print_pi_ij(int lx, int ly, double alpha, int i, int j){
    ewald_ND e2d(2,vector<int>({lx,ly}),alpha);
    double res = e2d.pi_ij(i,j);
    cout << res << '\n';
    return res;
    // cout << e2d.pi_1(i,j,e2d.kappa)+e2d.pi_2(i,j,e2d.kappa) << '\n';
}

// int main(){
//     int lx = 32, ly = 32, alpha = 3;
//     ewald_ND e2d(2,vector<int>({lx,ly}),alpha);
//     double J_tot = 0;
//     cout << "lx = 32, ly = 32 , alpha = 2" << '\n';
//     for(int j = 0; j < 32*32; j++){
//         for(int i = 0; i < 32*32; i++){
//             if(i == j) continue;
//             cout << j << " " << i << " " << e2d.pi_ij(i,j) << '\n';
//             J_tot += e2d.pi_ij(i,j);
//             // cout << j << " " << i << " " << e2d.dist_ij_min(j,i) << '\n';
//             // cout << e2d.xof(j) << " " << e2d.xof(i) << " " << '\n';
//             // cout << e2d.yof(j) << " " << e2d.yof(i) << " " << '\n';
//             // cout << j << " " << i << " " << pow(e2d.dist_ij_min(j,i),-alpha) << '\n';
//             // J_tot += pow(e2d.dist_ij_min(j,i),-alpha);
//         }
//     }
//     // cout << j << " " << i << " " << pow(e2d.dist_ij(i,j),-alpha) << '\n';
//     // int j = 0, i = 32;
//     // cout << j << " " << i << " ";
//     cout << J_tot << endl;
// }

// int main(){
//     int lx = 32, ly = 32, alpha = 3;
//     ewald_ND e2d(2,vector<int>({lx,ly}),alpha);
//     double J_tot = 0;
//     cout << "lx = 32, ly = 32 , alpha = 2" << '\n';
//     for(int j = 0; j < 32*32; j++){
//         for(int i = 0; i < 32*32; i++){
//             if(i == j) continue;
//             cout << j << " " << i << " " << e2d.pi_ij(i,j) << '\n';
//             J_tot += e2d.pi_ij(i,j);
//             // cout << j << " " << i << " " << e2d.dist_ij_min(j,i) << '\n';
//             // cout << e2d.xof(j) << " " << e2d.xof(i) << " " << '\n';
//             // cout << e2d.yof(j) << " " << e2d.yof(i) << " " << '\n';
//             // cout << j << " " << i << " " << pow(e2d.dist_ij_min(j,i),-alpha) << '\n';
//             // J_tot += pow(e2d.dist_ij_min(j,i),-alpha);
//         }
//     }
//     // cout << j << " " << i << " " << pow(e2d.dist_ij(i,j),-alpha) << '\n';
//     // int j = 0, i = 32;
//     // cout << j << " " << i << " ";
//     cout << J_tot << endl;
// }