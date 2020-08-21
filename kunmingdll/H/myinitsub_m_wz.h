#ifndef MYINITSUB_M_WZ_H
#define MYINITSUB_M_WZ_H

#include <cmath>

#include <QDebug>

#include "../writetofile.h"

void myInitSub_M_WZ(const double &newQ, const double &newD,
                    double **PipesectionTab, double **DenTab, double **VTab, int &FandDFlag,
                    int *FandDcount, const double *Pipelenth, const int &N, const double *ab_coef)
{
    double Diameter = 244.5;
    double Wallthickness = 9.49;
    double Diam = (Diameter - Wallthickness)/1000;
    double aup_dn1 = ab_coef[1];
    double bup_dn1 = ab_coef[2];
    double adn_up1 = ab_coef[3];
    double bdn_up1 = ab_coef[4];
    double aup_dn2 = ab_coef[5];
    double bup_dn2 = ab_coef[6];
    double adn_up2 = ab_coef[7];
    double bdn_up2 = ab_coef[8];
    double pi = 3.141592653589793;
    double u = newQ/(900*pi*pow(Diam,2));
    double newL = 60*u;
    int forepos1 = N - FandDcount[1] + 1;
    int flag_stop = 0;

    if(newQ == 0)
    {
        flag_stop = 1;
    }
    if(flag_stop == 0)
    {
        for(int i = forepos1; i<=N; i++)
        {
            DenTab[1][i-1] = DenTab[1][i];
            PipesectionTab[1][i-1] = PipesectionTab[1][i];

        }
    }
}

#endif // MYINITSUB_M_WZ_H
