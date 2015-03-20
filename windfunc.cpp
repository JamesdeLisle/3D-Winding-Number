#include "twodwind.h"
#define PI 3.141592

using namespace Eigen;

void chern::linspace(double pmin, double pmax, int pint)
{
    double step = (pmax-pmin) / (pint);
    while(pmin <= pmax)
    {
        pVec.push_back(pmin);
        pmin += step;
    }
}

chern::chern( double mu_, double delta_, int pint_ )
{
    int i,j,k;
    linspace(-PI,PI,pint_);
    double shift = 1e-2;

    //Matrix4cd Q, Qx, Qy, Qz, Id = Matrix4cd::Identity();
    //Matrix2cd q, qx, qy, qz;

    Matrix4cd Id = Matrix4cd::Identity();
    System.reserve(4);
    std::vector<Matrix2cd> q(4);
    std::vector<Matrix2cd> qpt(3);
    Matrix2cd qInv;
    double WindPart;
    std::vector<std::complex<double>> intpart(6);
    double gaptemp = 100.0;
    WindN = 0.0;

    for(i = 0; i < pint_; i++)
    {
        for( j = 0; j < pint_; j++ )
        {
            for ( k = 0; k < pint_; k++ )
            {


                System[0] = sys(pVec[i],pVec[j],pVec[k], mu_, delta_);
                System[1] = sys(pVec[i]+shift,pVec[j],pVec[k], mu_, delta_);
                System[2] = sys(pVec[i],pVec[j]+shift,pVec[k], mu_, delta_);
                System[3] = sys(pVec[i],pVec[j],pVec[k]+shift, mu_, delta_);
                for ( int l = 0; l < 4; l++)
                {
                    q[l] = (2.0 * (System[l].EVecOutput().col(0)*System[l].EVecOutput().col(0).adjoint() + System[l].EVecOutput().col(1)*System[l].EVecOutput().col(1).adjoint()) - Id).topRightCorner(2,2);
                }
                for ( int l = 0; l < 3; l++)
                {
                    qpt[l] = (q[1+l] - q[0])/shift;
                }
                qInv = q[0].inverse();
                intpart[0] = (qInv * qpt[0] * qInv * qpt[1] * qInv * qpt[2]).trace();
                intpart[1] = (qInv * qpt[1] * qInv * qpt[2] * qInv * qpt[0]).trace();
                intpart[2] = (qInv * qpt[2] * qInv * qpt[0] * qInv * qpt[1]).trace();
                intpart[3] = -(qInv * qpt[0] * qInv * qpt[2] * qInv * qpt[1]).trace();
                intpart[4] = -(qInv * qpt[2] * qInv * qpt[1] * qInv * qpt[0]).trace();
                intpart[5] = -(qInv * qpt[1] * qInv * qpt[0] * qInv * qpt[2]).trace();
                WindPart = std::real(intpart[0] + intpart[1] + intpart[2] + intpart[3] + intpart[4] + intpart[5]);
                WindN += WindPart;
                if ( System[0].getgap() < gaptemp )
                {
                    gaptemp = System[0].getgap();
                }
            }
        }
    }
    WindN = (1.0/(24.0*PI*PI)) * ((2.0 * PI)/pint_) * ((2.0 * PI)/pint_) * ((2.0 * PI)/pint_) * WindN;
    std::cout << WindN  << std::endl << std::endl;
    gapmin = gaptemp;
}


