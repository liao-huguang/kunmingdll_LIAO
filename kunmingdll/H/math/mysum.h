#ifndef MYSUM_H
#define MYSUM_H

double mySum(double*x, int N)
{
    double sum = 0;
    for(int i = 0; i<N; i++)
    {
        sum += x[i];
    }
    return sum;
}

#endif // MYSUM_H
