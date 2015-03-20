#include "twodwind.h"
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <eigen3/Eigen/Eigenvalues>

sys::sys(double pX_, double pY_, double pZ_, double mu_, double delta_)
{
    std::complex<double> ipz(0.0,pZ_);
    std::complex<double> e2ipz = exp(ipz);
    std::complex<double> e2mipz = exp(-ipz);
    std::complex<double> one(1.0,0.0);
    std::complex<double> I(0.0,1.0);
    std::complex<double> cosp( cos( (pX_ + pY_)/2.0 ), 0.0);
    std::complex<double> cosm( cos( (pX_ - pY_)/2.0 ), 0.0);

    double t = 4.0;
    double del = 2.0;
    std::complex<double> PHI = 2.0 * t *( e2ipz * cosp + e2mipz * cosm);
    std::complex<double> DELTA( 4.0 * del * ( cos( pX_/2.0 ) * cos( pY_/2.0 ) - cos( pY_/2.0 ) * cos( pZ_/2.0 ) - cos( pZ_/2.0 ) * cos( pX_/2.0 )), 0.0);
    std::complex<double> THETA( 4.0 * delta_ * cos( pZ_/2.0 ) * (cos( pY_/2.0 ) - cos( pX_/2.0 )) + mu_, 0.0);

    hamiltonian(0,0) = 0.0;
    hamiltonian(0,1) = 0.0;
    hamiltonian(0,2) = DELTA - I * THETA;
    hamiltonian(0,3) = -I * PHI;
    hamiltonian(1,1) = 0.0;
    hamiltonian(1,2) = -I * std::conj(PHI);
    hamiltonian(1,3) = DELTA + I * THETA;
    hamiltonian(2,2) = 0.0;
    hamiltonian(2,3) = 0.0;
    hamiltonian(3,3) = 0.0;

    hamiltonian(1,0) = std::conj(hamiltonian(0,1));
    hamiltonian(2,0) = std::conj(hamiltonian(0,2));
    hamiltonian(3,0) = std::conj(hamiltonian(0,3));
    hamiltonian(2,1) = std::conj(hamiltonian(1,2));
    hamiltonian(3,1) = std::conj(hamiltonian(1,3));
    hamiltonian(3,2) = std::conj(hamiltonian(2,3));

    hamiltonian = 0.5 * hamiltonian;

    Eigen::Matrix4cd HermT;

    HermT = hamiltonian - hamiltonian.adjoint();

    //std::cout << HermT;

    Eigen::ComplexEigenSolver<Eigen::Matrix4cd> ces(hamiltonian);

    Eigen::Vector4cd eigenvaltemp;
    eigenvaltemp = ces.eigenvalues();
    Eigen::Vector4d sorteigenvaltemp;

    Eigen::Vector4cd temp;

    sorteigenvaltemp[0] = std::real(eigenvaltemp[0]);
    sorteigenvaltemp[1] = std::real(eigenvaltemp[1]);
    sorteigenvaltemp[2] = std::real(eigenvaltemp[2]);
    sorteigenvaltemp[3] = std::real(eigenvaltemp[3]);
    //std::cout << sorteigenvaltemp;

    Eigen::Vector4d refer(0,1,2,3);
    int flag = 1;
    double tempeva, temprefer;
    for (int i = 0; (i < 4) && flag; i++)
    {
        flag = 0;
        for (int j = 0; j < 3; j++)
        {
            if (sorteigenvaltemp[j+1] < sorteigenvaltemp[j])
            {
                tempeva = sorteigenvaltemp[j];
                sorteigenvaltemp[j] = sorteigenvaltemp[j+1];
                sorteigenvaltemp[j+1] = tempeva;
                temprefer = refer[j];
                refer[j] = refer[j+1];
                refer[j+1] = temprefer;

                flag = 1;
            }
        }
    }


    for (int k = 0; k < 4; k++)
    {
        eigenvectors.col(k) << ces.eigenvectors().col(refer[k]);
    }

    eigenvalue = sorteigenvaltemp;
    gapdiff = std::real(eigenvalue[2] - eigenvalue[1]);
    //std::cout << gapdiff << " ";
}

