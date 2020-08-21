#ifndef MYQ_DTRANSSUUB_WZ_H
#define MYQ_DTRANSSUUB_WZ_H

#include <cmath>
#include <numeric>
#include <iostream>

#include "./math/mysum.h"

#include <QDebug>

void myQ_DTransSub_WZ(const double *fucQ, const double *fucD, const double *Q21, const double *D21, const int &Len, const int &t,
                      double **PipesectionTab, double **DenTab, double **VTab, int &FandDFlag, int *FandDCount, const double *Pipelength,
                      const int &N, const double *ab_coef)
{
    double pi = 3.141592653589793;

    double *Q = new double[Len+1];
    doubel *D = new double[Len+1];
    for(int i=1; i<=Len; i++)
    {
        Q[i] = fucQ[i];
        D[i] = fucD[i];
    }

    double Diameter = 244.5;
    double Wallthickness = 9.79;
    double Diam = (Diameter - Wallthickness - Wallthickness)/1000;
    double aup_dn1 = ab_coef[1];
    double bup_dn1 = ab_coef[2];
    double adn_up1 = ab_coef[3];
    double bdn_up1 = ab_coef[4];
    double aup_dn2 = ab_coef[5];
    double bup_dn2 = ab_coef[6];
    double adn_up2 = ab_coef[7];
    double bdn_up2 = ab_coef[8];

    double *flow = new double[1440+1+t];
    double *conc = new double[1440+1+t];
    for(int i = 1; i <= 1440; i++)
    {
        flow[i] = Q[i];
        conc[i] = D[i];
    }
    for(int i = 1441; i<=1440+t; i++)
    {
        flow[i] = Q21[i-1440];
        conc[i] = D21[i-1440];
    }
    int time = t + Len;

    double *len_L1 = new double[Len+1];
    double *len_L2 = new double[Len+1];
    double *conc_L1 = new double[Len+1];
    double *conc_L2 = new double[Len+1];
    double *vup_L1 = new double[Len+1];
    double *vdn_L1 = new double[Len+1];
    double *vup_L2 = new double[Len+1];
    double *vdn_L2 = new double[Len+1];
    for(int i = 1; i<=Len; i++)
    {
        len_L1[i] = 0;
        len_L2[i] = 0;
        conc_L1[i] = 0;
        conc_L2[i] = 0;
        vup_L1[i] = 0;
        vdn_L1[i] = 0;
        vup_L2[i] = 0;
        vdn_L2[i] = 0;
    }
    int temp1 = 1;
    int temp2 = 1;
    for(int i=time; i>=1; i--)
    {
        if(mySum(len_L1+1, Len) < Pipelength[1])
        {
            if (flow[i] == 0)
            {
                len_L1[temp1] = 100;
            }
            else
            {
                len_L1[temp1] = flow[i]/(15*pi*Diam*Diam);
            }
            conc_L1[temp1] = conc[i];
            if (mySum(len_L1+1, Len) >= Pipelength[1])
            {
                double temp_sum = mySum(len_L1+1, Len);
                len_L1[temp1] = len_L1[temp1] - (temp_sum - Pipelength[1]);
                len_L2[temp2] = temp_sum - Pipelength[1];
                conc_L2[temp2] = conc[i];
                temp2 = temp2 + 1;
            }
            temp1 = temp1 + 1;
        }
        else
        {
            if (flow[i] == 0)
            {
                len_L2[temp2] = 100;
            }
            else
            {
                len_L2[temp2] = flow[i];
            }
            conc_L2[temp2] = conc[i];
            temp2 = temp2 + 1;
            if (mySum(len_L2+1, Len) >= Pipelength[2])
            {
                double temp_sum = mySum(len_L2+1, Len);
                len_L2[temp2-1] = len_L2[temp2-1] - (temp_sum - Pipelength[2]);
                break;
            }
        }
    }
    FandDCount[1] = temp1 - 1;
    FandDcount[2] = temp2 - 1;

    for(int i = 1; i <= Len; i++)
    {
        vup_L1[i] = sqrt(aup_dn1*(bup_dn1 - conc_L1[i]));
        vdn_L1[i] = sqrt(adn_up1*(bdn_up1 - conc_L1[i]));
        vup_L2[i] = sqrt(aup_dn2*(bup_dn2 - conc_L2[i]));
        vdn_L2[i] = sqrt(adn_up2*(bdn_up2 - conc_L2[i]));
    }
    for(int i = 1; i<=FandDCount[1]; i++)
    {
        PipesectionTab[1][N-i+1] = len_L1[i];
        DenTab[1][N-i+1] = conc_L1[i];
        VTab[1][N-i+1] = vup_L1[i];
        VTab[2][N-i+1] = vdn_L1[i];
    }
    for(int i=1; i<=FandDCount[2]; i++)
    {
        PipesectionTab[2][N-i+1] = len_L2[i];
        DenTab[2][N-i+1] = conc_L2[i];
        VTab[3][N-i+1] = vup_L2[i];
        VTab[4][N-i+1] = vdn_L2[i];
    }

    FandDFlag = 1;

    delete [] len_L1;
    delete [] len_L2;
    delete [] conc_L1;
    delete [] conc_L2;
    delete [] vup_L1;
    delete [] vdn_L1;
    delete [] vdn_L2;
    delete [] vup_L2;

    delete [] flow;
    delete [] conc;
    delete [] Q;
    delete [] D;
}

#endif // MYQ_DTRANSSUUB_WZ_H
