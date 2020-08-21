#ifndef LEAKDIAGNOSE_H
#define LEAKDIAGNOSE_H

#include <QtCore/qglobal.h>

#include <QString>

__declspec(dllexport) int leakDiagnose(const bool &fidl, const bool &fid2, double **ZQ, const int &NL, const double *density,
                                       int &TFlagP1, double **P1, double **SP1, int &OldCnt1, double **OldAlarmTab1, const int &TabLen,
                                       int &TFlagP2, double **P2, double **SP2, int &OldCnt2, double **OldAlarmTab2,
                                       const double *Q1, const double *D1, const int &timeCount,
                                       double **PipesectionTab, double **DenTab, double **VTab, int &FandDFlag, int *FandDcount,
                                       const bool &SingleSignalLeakDetectionFlag, const bool &stationAlarmFlag,
                                       QString *AlarmLevel, double *leakLocation, int *signalPosx, int *signalPosy, double *myPeak,
                                       bool *pipeLineDiagnoseOpenFlag);

#endif // LEAKDIAGNOSE_H
