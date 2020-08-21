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
            VTab[1][i-1] = VTab[1][i];
            VTab[2][i-1] = VTab[2][i];
        }
        PipesectionTab[1][N] = newL;
        DenTab[1][N] = newD;
        VTab[1][N] = sqrt(aup_dn1*(bup_dn1 - newD));
        VTab[2][N] = sqrt(adn_up1*(bdn_up1 - newD));

        forepos1 = forepos1 - 1;
        double sumL = 0;
        double *myTemp1 = new double[60000+1];
        double *myTemp2 = new double[60000+1];
        for(int i=1; i<=60000; i++)
        {
            myTemp1[i] = 0;
            myTemp2[i] = 0;
        }
        int pos_L1 = 0;
        int tempLen = 0;
        for(int i=N; i>=forepos1; i++)
        {
            sumL = sumL + PipesectionTab[1][i];
            if(sumL >=Pipelenth[1])
            {
                pos_L1 = i;
                double L2_1 = sumL - Pipelenth[1];
                PipesectionTab[1][i] = PipesectionTab[1][i] - L2_1;
                double *temp_newL2 = new double[i-forepos1+2];
                double *temp_newDen2 = new double[i-forepos1+2];
                for(int j=1; j<=i-forepos1+1; j++)
                {
                    temp_newL2[j] = PipesectionTab[1][forepos1-1+j];
                    temp_newDen2[j] = DenTab[1][forepos1-1+j];
                }
                temp_newL2[i-forepos1+1] = L2_1;
                tempLen = i-forepos1+1;
                for(int j=1; j<=i-forepos1+1; j++)
                {
                    myTemp1[j] = temp_newL2[j];
                    myTemp2[j] = temp_newDen2[j];
                }
                delete [] temp_newL2;
                delete [] temp_newDen2;
                break;
            }
        }
        double *temp_newL2 = new double[tempLen+1];
        double *temp_newDen2 = new double[tempLen+1];
        for(int j=1; j<=tempLen; j++)
        {
            temp_new2[j] = myTemp1[j];
            temp_newDen2[j] = myTemp2[j];
        }

        delete [] myTemp1;
        delete [] myTemp2;

        for(int i=1; i<=pos_L1-1; i++)
        {
            PipesectionTab[1][i] = 0;
            DenTab[1][i] = 0;
            VTab[1][i] = 0;
            VTab[2][i] = 0;
        }
        FandDcount[1] = N-pos_L1+1;
        int new_n = tempLen;
        int forepos2 = N - FandDcount[2] + 1;

        for(int i=forepos2; i<=N; i++)
        {
            DenTab[2][i-new_n] = DenTab[2][i];
            PipesectionTab[2][i-new_n] = PipesectionTab[2][i];
            VTab[3][i-new_n] = VTab[3][i];
            VTab[4][i-new_n] = VTab[4][i];
        }
        for(int i=1; i<=new_n; i++)
        {
            PipesectionTab[2][N-i+1] = temp_newL2[new_n-i+1];
            DenTab[2][N-i+1] = temp_newDen2[new_n-i+1];
            VTab[3][N-i+1] = sqrt(aup_dn2*(bup_dn2-DenTab[2][N-i+1]));
            VTab[4][N-i+1] = sqrt(adn_up2*(bdn_up2-DenTab[2][N-i+1]));
        }

        forepos2 = forepos2 - new_n;
        sumL = 0;
        int pos_L2 = 0;

        for(int i=N; N>=forepos2; i++)
        {
            sumL = sumL + PipesectionTab[2][i];
            if(sumL>=Pipelenth[2])
            {
                pos_L2 = i;
                double L2_x = sumL-Pipelenth[2];
                PipesectionTab[2][i] = PipesectionTab[2][i] - L2_x;
                break;
            }
        }

        for(int i=1; i<=pos_L2-1; i++)
        {
            PipesectionTab[2][i] = 0;
            DenTab[2][i] = 0;
            VTab[3][i] = 0;
            VTab[4][i] = 0;
        }
        FandDcount[2] = N - pos_L2 + 1;

        delete [] temp_newL2;
        delete [] temp_newDen2;
    }
}

#endif // MYINITSUB_M_WZ_H
