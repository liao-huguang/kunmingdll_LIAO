#include "kunmingdll.h"

#include <QString>

#include <cmath>

#include "./H/myptransfunc_lp.h"
#include "./H/mymainfunc.h"
#include "./H/myinitsub_m.h"

#include "./H/myinitsub_m_wz.h"
#include "./H/myq_dtranssub.h"
#include "./H/myq_dtranssub_wz.h"

int leakDiagnose(const bool &fidl, const bool &fid2, double **ZQ, const int &NL, const double *density,
                 int &TFlagP1, double **P1, double **SP1, int &OldCnt1, double **OldAlarmTab1, const int &TabLen,
                 int &TFlagP2, double **P2, double **SP2, int &OldCnt2, double **OldAlarmTab2,
                 const double *Q1, const double *D1, const int &timeCount,
                 double **PipesectionTab, double **DenTab, double **VTab, int &FandDFlag, int *FandDcount,
                 const bool &SingleSignalLeakDetectionFlag, const bool &stationAlarmFlag,
                 QString *AlarmLevel, double *leakLocation, int *signalPosx, int *signalPosy, double *myPeak,
                 bool *pipeLineDiagnoseOpenFlag)
{
    AlarmLevel[1].clear();
    AlarmLevel[2].clear();
    for(int i = 1; i<=2; i++)
    {
        AlarmLevel[i].clear();
        leakLocation[i] = -99;
        signalPosx[i] = -99;
        signalPosy[i] = -99;
        myPeak[i] = 0;
    }
    double Params11[16+1] = {0, 2000*pow(10,-12), 2000*pow(10,6), 20000.0, 0.75, 7, 20000, 0.75, 1.2, 4, 200, 0.05, 80, 0.001, 20, 0.001, 750};
    double Params12[11+1] = {0, 47.410, 100.0, 3.5, -3.5, 3.5, -3.5, 1108.6, 1080.5, 5000, 1300, 1080.5};
    double Params21[16+1] = {0, 2000*pow(10,-12), 2000*pow(10,6), 47000.0, 0.75, 7, 20000, 0.75, 1.2, 4, 200, 0.05, 40, 0.001, 10, 0.001, 750};
    double Params22[11+1] = {0, 41.834, 100.0, 3.5, -3.5, 3.5, -3.5, 987.2563, 987.2563, 4700, 1245, 987.2563};

    double Pipelength[2+1] = {0, Params12[1]*1000, Params22[1]*1000};

    const int Nd = 1440;
    const int QCnt = Nd + timeCount;

    const int N = 6000;
    const int MLen = 10;
    const int Dlen = MLen*N;
    const int CH = 2;

    double ab_coef[8+1] = {0, 14007.3710120125, 152.274296132468, 13025.0686600743, 155.794107925798, 10555.4116138042, 162.692190990793, 10555.4116138042, 162.692190990793};

    double *Q = new double[Nd+1];
    double *D = new double[Nd+1];
    double *Q21 = new double[Nd+1];
    for(int i = 1; i<=Nd; i++)
    {
        Q[i] = Q1[i];
        D[i] = D1[i];

        Q21[i] = 0;
        D21[i] = 0;
    }
    for(int i=1; i<=timeCount; i++)
    {
        Q21[i] = Q1[i+Nd];
        D21[i] = D1[i+Nd];
    }

    if(FandDFlag == 0)
    {
        myQ_DTransSub_WZ(Q,D,Q21,D21,Nd,timeCount,PipesectionTab, DenTab, VTab, FandDFlag, FandDcount, Pipelength, DLen, ab_coef);
    }
    else
    {
        myInitSub_M_WZ(Q21[timeCount],D21[timeCount], PipesectionTab, DenTab, VTab, FandDcount, Pipelength, Dlen, ab_coef);
    }


    if(pipeLineDiagnoseOpenFlag[1])
    {
        if(fid1)
        {

        }
        else
        {
            for(int j = 1; j <= N; j++)
            {
                P1[1][j] = 0;
                p1[2][j] = 0;
            }
            TFlagP1 = 0;
        }
        if(TFlagP1>=MLen)
        {
            myMainFunc(TabLen, OldAlarmTab1, OldCnt1, ZQ, SP1, CH, DLen, N, Params12, stationAlarmFlag, SingleSignalLeakDetectionFlag,
                       PipesectionTab[1], DenTab[1], const_cast<const double**>(VTab), FandDcount[1],
                       AlarmLevel[1], leakLocation[1], signalPosx[1], signalPosy[1], myPeak[1], 1);
        }
    }

    if(pipeLineDiagnoseOpenFlag[2])
    {
        if(fid2)
        {

        }
        else
        {
            for(int j = 1; j <= N; j++)
            {
                P2[1][j] = 0;
                P2[2][j] = 0;
            }
            TFlagP2 = 0;
        }
        if(TFlagP2>=MLen)
        {
            myMainFunc(TabLen, OldAlarmTab2, OldCnt2, ZQ, SP2, CH, DLen, N, Params21, Params22, stationAlarmFlag, SingleSignalLeakDetectionFlag,
                       PipesectionTab[2], DenTab[2], const_cast<const double**>(VTab+2), FandDcount[2],
                       AlarmLevel[2], leakLocation[2], signalPosx[2], signalPosy[2], myPeak[2], 2);
        }
    }
    delete [] Q;
    delete [] D;
    delete [] Q21;
    delete [] D21;
    return 0;
}
