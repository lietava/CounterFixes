#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h> 
#include <TCanvas.h>
#include <TGraph.h>
#include <TApplication.h>
#include <cmath>
#include <vector>
#include <iostream>

const double alpha = 2e-6;
const double L = 2.0;
const int N = 100;
const double dx = L / (N - 1);
const double dt = 100;
const double r = alpha * dt / (dx * dx);
const int steps = 3600;
double LB = 10;
double RB = 20;

double fun(double x)
{
    //return (10. + std::sin(M_PI * x));
    return 5*x + 10;
    //return x+10;
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    TCanvas* c1 = new TCanvas("c1", "Crank–Nicolson Result", 800, 600);
    //TGraph* g = new TGraph(N-2);
    std::array<TGraph*,steps+1> tgraphs;
    float *xx = new float[steps];
    float *yy = new float[steps]; 
    float *yyG = new float[steps];
    // Initialize u(x,0) = sin(pi*x)
    TVectorD u(N-2), u_new(N-2);
    for (int i = 0; i < N-2; ++i) {
        double x = (i+1) * dx;
        u[i] = fun(x);
    }
    LB=10;
    RB=20;
    // Tridiagonal matrix A and RHS vector b
    TMatrixD A(N-2, N-2);
    TVectorD b(N-2);

    // Fill matrix A (Crank–Nicolson implicit matrix)
    A(0,0) = (2*r+1);
    A(0,1) =  -r;
    for (int i = 1; i < N-3; ++i) {
        A(i, i)   = 2*r+1;
        A(i, i-1) = -r;
        A(i, i+1) = -r;
    }
    A(N-3,N-3) = 2*r+1;
    A(N-3,N-4) = -r;
    //
    tgraphs[0] = new TGraph(N-2);
    for (int i = 0; i < N-2; ++i) {
        tgraphs[0]->SetPoint(i, (i+1) * dx, u[i] );
    }
    // Time stepping
    for (int t = 0; t < steps; ++t) {
        // Build RHS vector b
        RB = 20+10*sin(t*dt*2*3.14/3600./24.);
        for (int i = 1; i < N-3; ++i) {
            b[i] = u[i];
        }
        b[0] = u[0] + r * LB;
        b[N-3] = u[N-3] + r * RB;
        // Solve the linear system A * u_new_inner = b
        TDecompLU solver(A);
        //solver.Print();
        //std::cout << "RHS b:" << std::endl;
        //b.Print();
        bool ok;
        TVectorD u_new = solver.Solve(b,ok);
        if(!ok) {
            std::cout << "===> NOT OK" << std::endl;
            return 1;
        }

        // Set boundary conditions and update full vector
        u = u_new;
        //u.Print();
        
        tgraphs[t+1] = new TGraph(N-2);
        for (int i = 0; i < N-2; ++i) {
            tgraphs[t+1]->SetPoint(i, (i+1) * dx, u[i] );
        }
        xx[t] = t*dt;
        yy[t] = u[N*1.2/2.];
        yyG[t] = RB;
    }

    // Plot the final result
    TGraph* bcon = new TGraph[2];
    bcon->SetPoint(0,0,LB);
    bcon->SetPoint(1,L,RB);
    bcon->SetTitle("Heat Equation Solution;Position x;Temperature u(x)");
    bcon->SetMarkerStyle(21);
    //bcon->Draw("AP");
    for(int i = 0;  i < steps+1; i++) {
        tgraphs[i] ->Draw("L SAME");
    }
    double chi = 0;
    double x = 0;
    for(int i = 0; i < N-2; i++) {
        chi += (30*(i+1)*dx - u[i])*(30*(i+1)*dx - u[i]);
    }
    std::cout << "chi:" << chi << std::endl;
    TGraph* g1 = new TGraph(steps,xx,yy);
    g1->SetMarkerStyle(21);
    TGraph* g2 = new TGraph(steps,xx,yyG);
    g2->Draw("AP");
    g1->Draw("SAME");
    app.Run();
    return 0;
}