///////Devendra Chaudhari 
//////couette Flow
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
// relaxation time calculation
double Relaxation_time(double Re, double Ma, int L)
{
    double U, kin_visc, relax_time;
    double cs = 1 / sqrt(3.0);
    U = (Ma * cs);
    kin_visc = U * L / Re;
    relax_time = 3 * kin_visc + 0.5;
    return relax_time;
}

int main()
{
    double Re = 100, Ma = 0.075, tau, ux, uy, rho, u, term1, term2;
    int Nx, Ny, ix, iy;
    const int q = 9;
    const double rho0 = 1.0;
    double cs = 1 / sqrt(3.0);
    int iteration = 0;
    double error = 0;
    double uw = Ma*cs;
    std::cout << "enter number of nodes in y-direction: ";
    std::cin >> Ny;
    Nx= Ny;

    std::vector<std::vector<double>> ux_old(Nx + 3, std::vector<double>(Ny + 3, 1.0));
    std::vector<std::vector<double>> uy_old(Nx + 3, std::vector<double>(Ny + 3, 1.0));
    std::vector<std::vector<double>> ux_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> uy_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> u_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> rho_final(Nx + 3, std::vector<double>(Ny + 3, rho0));
    std::vector<std::vector<double>> p(Nx + 3, std::vector<double>(Ny + 3, 0.0));

    std::setprecision(7);

    tau = Relaxation_time(Re, Ma, (Ny));

    int cx[q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    double w[q] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
    int opposite[q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    std::vector<std::vector<std::vector<double>>> f(q, std::vector<std::vector<double>>(Nx + 3, std::vector<double>(Ny + 3, 0.0)));
    std::vector<std::vector<std::vector<double>>> fcoll(q, std::vector<std::vector<double>>(Nx + 3, std::vector<double>(Ny + 3, 0.0)));
    std::vector<double> feq(q, 0.0);

    // Initialization polpulation

    for (int i = 1; i <= Nx + 1; i++)
    {
        for (int j = 1; j <= Ny + 1; j++)
        {
            for (int a = 0; a < q; ++a)
            {
                f[a][i][j] = w[a] * rho0;
            }
        }
    }

    std::ofstream file;
    file.open("lbm_output.dat");

    std::ofstream data;
    data.open("error.txt");
    data << "x-grid ponts" << Nx + 1 << "\t" << "Y-grid points" << Ny + 1 << "\n";
    data << "i" << "j" << "\t\t\t" << "rho\t\t\t" << "ux\t\t\t" << "uy\t\t\t" << "u\t\t\t" << "psi\t\t\t" << "omega\t\t\t" << "error\n";
   

    std::ofstream ux_vel_5;
    ux_vel_5.open("ux_velocity5.dat");

    do
    {
        //////////////////collision///////////////////////////////////////
        for (int i = 1; i <= Nx + 1; i++)
        {
            for (int j = 1; j <= Ny + 1; j++)
            {
                ux = ux_final[i][j];
                uy = uy_final[i][j];
                rho = rho_final[i][j];
                for (int a = 0; a < q; a++)
                {
                    term1 = ux * cx[a] + uy * cy[a];
                    term2 = term1 * term1;
                    feq[a] = w[a] * rho * (1.0 + 3.0 * term1 + 4.5 * term2 - 1.5 * (ux * ux + uy * uy));
                    fcoll[a][i][j] = f[a][i][j] - (f[a][i][j] - feq[a]) / tau; // collision
                }
            }
        }

        ////////////////////////streaming/////////////////////////
        for (int i = 1; i <= Nx + 1; i++)
        {
            for (int j = 1; j <= Ny + 1; j++)
            {
                for (int a = 0; a < q; a++)
                {
                    int nextX = i + cx[a];
                    int nextY = j + cy[a];
                    f[a][nextX][nextY] = fcoll[a][i][j];
                }
            }
        }

        ///////Boundary Conditions
        double uxt = uw, uyt = 0.0, uxb = 0.0, uyb = 0.0, uxl=0.0 , uyl = 0.0, uyr=0.0,uxr=0.0;

        // Top Wall  ux=uw uy=0.0 rho= unknown.
        for (int i = 1; i <= Nx + 1; i++)
        {
            
            rho_final[i][Ny + 1] = 1.0;
            f[4][i][Ny + 1] = f[2][i][Ny + 1] - ((2.0 / 3.0) * (rho_final[i][Ny + 1]) * (uyt));
            f[7][i][Ny + 1] = f[5][i][Ny + 1] + (0.5 * (f[1][i][Ny + 1] - f[3][i][Ny + 1])) - ((1.0 / 6.0) * (rho_final[i][Ny + 1] * uyt)) - (0.5 * (rho_final[i][Ny + 1] * uxt));
            f[8][i][Ny + 1] = f[6][i][Ny + 1] + (0.5 * (f[3][i][Ny + 1] - f[1][i][Ny + 1])) - ((1.0 / 6.0) * (rho_final[i][Ny + 1] * uyt)) + (0.5 * (rho_final[i][Ny + 1] * uxt));
        }

        // Bottom Wall ux=0 uy=0 rho= unknown
        for (int i = 1; i <= Nx + 1; i++)
        {
           
            rho_final[i][1]=1.0;
            f[2][i][1] = f[4][i][1] + ((2.0 / 3.0) * rho_final[i][1] * uyb);
            f[5][i][1] = f[7][i][1] - (0.5 * (f[1][i][1] - f[3][i][1])) + ((1.0 / 6.0) * (rho_final[i][1] * uyb)) + (0.5 * (rho_final[i][1] * uxb));
            f[6][i][1] = f[8][i][1] - (0.5 * (f[3][i][1] - f[1][i][1])) + ((1.0 / 6.0) * (rho_final[i][1] * uyb)) - (0.5 * (rho_final[i][1] * uxb));
        }

        ////periodic Bc
        //left wall
        for(int j=2;j<=Ny;j++)
        {
            f[1][1][j] = f[1][Nx+1][j] ;
            f[5][1][j] = f[5][Nx+1][j] ;
            f[8][1][j] = f[8][Nx+1][j] ;
        }

        ///Right Wall
        for(int j=2;j<=Ny;j++)
        {
            f[3][Nx + 1][j] = f[3][1][j] ;
            f[7][Nx + 1][j] = f[7][1][j] ;
            f[6][Nx + 1][j] = f[6][1][j] ;
        }
 
        // Corner treatment

        // bottom left
        // concave
        f[1][1][1] = f[3][1][1] + ((2.0 / 3.0) * rho_final[1][2] * uxb);
        f[2][1][1] = f[4][1][1] + ((2.0 / 3.0) * rho_final[1][2] * uxb);
        f[5][1][1] = f[7][1][1] + ((1.0 / 6.0) * rho_final[1][2] * (uxb + uyb));
        f[6][1][1]=((1.0/12.0)*rho_final[1][2]*(uyb-uxb));
        f[8][1][1]=((1.0/12.0)*rho_final[1][2]*(uxb-uyb));
        f[0][1][1]= rho_final[1][2]- (f[1][1][1]+ f[2][1][1]+f[3][1][1]+f[4][1][1]+f[5][1][1]+f[6][1][1]+f[7][1][1]+f[8][1][1]);
        
        // Top Left

        f[1][1][Ny + 1] = f[3][1][Ny + 1] + ((2.0 / 3.0) * rho_final[1][Ny] * uxt);
        f[4][1][Ny + 1] = f[2][1][Ny + 1] - ((2.0 / 3.0) * rho_final[1][Ny] * uyt);
        f[8][1][Ny + 1] = f[6][1][Ny + 1] + ((1.0 / 6.0) * rho_final[1][Ny] * (uxt - uyl));
        f[5][1][Ny+1]=((1.0/12.0)*rho_final[1][Ny]*(uxt+uyt));
        f[7][1][Ny+1]=((1.0/12.0)*rho_final[1][Ny]*(-uxt-uyt));
        f[0][1][Ny+1]= rho_final[1][Ny]- (f[1][1][Ny+1]+ f[2][1][Ny+1]+f[3][1][Ny+1]+f[4][1][Ny+1]+f[5][1][Ny+1]+f[6][1][Ny+1]+f[7][1][Ny+1]+f[8][1][Ny+1]);
      
        // top Right
        
        f[3][Nx + 1][Ny + 1] = f[1][Nx + 1][Ny + 1] - ((2.0 / 3.0) * rho_final[Nx + 1][Ny] * uxt);
        f[4][Nx + 1][Ny + 1] = f[2][Nx + 1][Ny + 1] - ((2.0 / 3.0) * rho_final[Nx + 1][Ny] * uyt);
        f[7][Nx + 1][Ny + 1] = f[5][Nx + 1][Ny + 1] - ((1.0 / 6.0) * rho_final[Nx + 1][Ny] * (uxt + uyt));
        f[6][Nx+1][Ny+1]=((1.0/12.0)*rho_final[Nx+1][Ny]*(uyt-uxt));
        f[8][Nx+1][Ny+1]=((1.0/12.0)*rho_final[Nx+1][Ny]*(uxt-uyt));
        f[0][Nx+1][Ny+1]= rho_final[Nx+1][Ny]- (f[1][Nx+1][Ny+1]+ f[2][Nx+1][Ny+1]+f[3][Nx+1][Ny+1]+f[4][Nx+1][Ny+1]+f[5][Nx+1][Ny+1]+f[6][Nx+1][Ny+1]+f[7][Nx+1][Ny+1]+f[8][Nx+1][Ny+1]);
       

        // Bottom right
        // concave
        f[3][Nx + 1][1] = f[1][Nx + 1][1] - ((2.0 / 3.0) * rho_final[Nx + 1][2] * uxb);
        f[2][Nx + 1][1] = f[4][Nx + 1][1] + ((2.0 / 3.0) * rho_final[Nx + 1][2] * uyb);
        f[6][Nx + 1][1] = f[8][Nx + 1][1] - ((1.0 / 6.0) * rho_final[Nx + 1][2] * (uxb - uyb));
        f[5][Nx+1][1]=((1.0/12.0)*rho_final[Nx+1][2]*(uxb+uyb));
        f[7][Nx+1][1]=((1.0/12.0)*rho_final[Nx+1][2]*(-uxb-uyb));
        f[0][Nx+1][1]= rho_final[Nx+1][2]- (f[1][Nx+1][1]+ f[2][Nx+1][1]+f[3][Nx+1][1]+f[4][Nx+1][1]+f[5][Nx+1][1]+f[6][Nx+1][1]+f[7][Nx+1][1]+f[8][Nx+1][1]);
     
    
        /////////////////////macros////////////////
        double sum_unew = 0.;
        double sum_uold = 0.;
        for (int i = 1; i <= Nx + 1; i++)
        {
            for (int j = 1; j <= Ny + 1; j++)
            {
                ux = 0.0;
                uy = 0.0;
                rho = 0.0;
                for (int a = 0; a < q; a++)
                {
                    rho += f[a][i][j];
                    ux += f[a][i][j] * cx[a];
                    uy += f[a][i][j] * cy[a];
                }
                ux /= rho;
                uy /= rho;
                ux_final[i][j] = ux;
                uy_final[i][j] = uy;
                double t = (ux * ux) + (uy * uy);
                u_final[i][j] = sqrt(t);
                rho_final[i][j] = rho;
                p[i][j] = ((rho) * (u_final[i][j]));

                sum_unew += (ux - ux_old[i][j]) * (ux - ux_old[i][j]) + (uy - uy_old[i][j]) * (uy - uy_old[i][j]);
                sum_uold += ux_old[i][j] * ux_old[i][j] + uy_old[i][j] * uy_old[i][j];

                ux_old[i][j] = ux;
                uy_old[i][j] = uy;
            }
        }
        error = sqrt(sum_unew / sum_uold);

        if (iteration % 100 == 1)
        {
            std::cout << "Iteration " << iteration << ": Error = " << std::setprecision(7) << error << "\n";
        }

        iteration++;

    } while (error > 1e-06);
    std::vector<std::vector<double>> psi(Nx + 3, std::vector<double>(Ny + 3, 0.0));


    for (int i = 1; i <= Nx+1; i++) {
        psi[i][1] = psi[i-1][1] - uy_final[i][1] * 1.0;
    }

    // Then, integrate along the y-direction for each column
    for (int i = 1; i < Nx+1; i++) 
    {
        for (int j = 2; j < Ny+1; ++j) 
        {
            psi[i][j] = psi[i][j-1] + ux_final[i][j] * 1.0;
        }
    }

    // std::vector<std::vector<double>> omega(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    // // Calculate vorticity
    // for (int i = 1; i < Nx+1; i++) {
    //     for (int j = 1; j < Ny+1; j++) {
    //         // Calculate the Laplacian of psi using central differences
    //         double d2psi_dx2 = (psi[i+1][j] - 2 * psi[i][j] + psi[i-1][j]) ;
    //         double d2psi_dy2 = (psi[i][j+1] - 2 * psi[i][j] + psi[i][j-1]) ;

    //         // Calculate vorticity
    //         omega[i][j] = -(d2psi_dx2 + d2psi_dy2);
    //     }
    // }
    for (int i = 1; i <= Nx + 1; i++)
    {
        for (int j = 0; j <= Ny + 1; j++)
        {
            data << i << j << "\t\t" << rho_final[i][j] << "\t\t" << ux_final[i][j] << "\t\t" << uy_final[i][j] << "\t\t" << u_final[i][j] << "\t\t" << psi[i][j] << "\t\t"  << "\n";
        }
    }

    file << "TITLE = \"couette\"\n";
    file << "VARIABLES = \"X\", \"Y\", \"Ux\", \"Uy\", \"Rho\",\"psi\", \"p\"\n";
    file << "ZONE I=" << Nx + 1 << ", J=" << Ny + 1 << ", F=POINT\n";

    for (int i = 1; i <= Nx + 1; i++)
    {
        for (int j = 1; j <= Ny + 1; j++)
        {
            file << i << " " << j << " " << ux_final[i][j] << " " << uy_final[i][j] << " " << rho_final[i][j] * (1. / 3.) << " " << psi[i][j] << " " << " " << p[i][j] << " " << "\n";
        }
    }
    
    for (int j = 1; j <= Ny + 1; j++)
    {
        double i2 = (double)Ny;
        double yaxis = (j) / i2;
        
        ux_vel_5 << "\t" << (ux_final[((Nx-1)/2.0)][j] / uw) << "\t" << yaxis << "\n";
    }

    std::vector<double>ux_exact(Nx+3,0.0);
    
    std::ofstream exact;
    exact.open("ux_exact.dat");

    for(int i=1;i<=Ny+1;i++)
    {
        double i3 = (double)Ny;
        double yaxis1 = (i) / i3;
        ux_exact[i]= uw*((double) (i-1)/(double) Ny);
        exact<<(ux_exact[i]/uw)<<"\t"<<yaxis1<<"\n";
    }
    

    ux_vel_5.close();
    exact.close();
    file.close();
    data.close();
    
    return 0;
}
