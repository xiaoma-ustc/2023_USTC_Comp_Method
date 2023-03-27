#include "func.h"
#include <math.h>
#include <iomanip>
Solution::Solution(int n)
{
    A.resize(n, std::vector<double>(n));
    zero.resize(n);
    b.resize(n);
    y.resize(n);
}

void Solution::GaussEliminationWithPivoting(double epsilon)
{
    InitVector(0.01*0.01*0.5);
    b[99] -= (epsilon + 0.01);
    int n = A.size();
    for(int i = 0; i < n - 1; ++i)
    {
        if(abs(A[i][i]) < abs(A[i + 1][i]))
        {
            for(int j = i; j < n; ++j)
            {
                std::swap(A[i][j], A[i + 1][j]);
                std::swap(b[i], b[i + 1]);
            }
        }
        double t= A[i + 1][i] / A[i][i];
        A[i + 1][i] = 0;
        for(int j = i + 1; j < n; ++j)
        {
            A[i + 1][j] -= t * A[i][j];
        }
        b[i + 1] -= t * b[i];
    }
    b[n - 1] /= A[n - 1][n - 1];
    for(int i = n - 2; i >= 0; --i)
    {
        for(int j = i + 1; j < n; ++j)
        {
            b[i] -= b[j] * A[i][j];
        }
        b[i] /= A[i][i];
    }

    std::cout<<"Result of Guass Elimination With Pivoting "<<std::endl;
    for(int i = 0; i < n; ++i)
    {
        std::cout<<std::setw(10)<<b[i]<<"\t";
        if(i > 1 && (i + 1) % 5 == 0)
        {
            std::cout<<std::endl;
        }
    }
    std::cout<<std::endl;
    std::cout<<std::left<<"Acc : "<<Distance(b, y)/Distance(y, zero)<<std::endl;
    b = std::vector<double>(b.size());

    
}

void Solution::GaussSeidelIteration(double epsilon)
{
    int n = A.size();
    InitVector(0.01*0.01*0.5);
    b[99] -= (epsilon + 0.01);
    std::vector<double> x1(n);
    std::vector<double> x2 = x1;
    double e = 0.000001;
    do
    {
        for(int i = 1; i < n; ++i)
        {
            x1 = x2;
            for(int j = 0; j < n; ++j)
            {
                double s = 0;
                for(int k = 0; k < n; ++k)
                {
                    s += A[j][k] * x2[k];
                }
                x2[j] = (b[j] - s + A[j][j] * x2[j]) / A[j][j];
            }
        }
    } while (JudgeIteration(x1, x2, e));
    
    std::cout<<"Result of Guass Seidel Iteration "<<std::endl;
    for(int i = 0; i < n; ++i)
    {
        std::cout<<std::setw(10)<<x2[i]<<"\t";
        if(i > 1 && (i + 1) % 5 == 0)
        {
            std::cout<<std::endl;
        }
    }
    std::cout<<std::endl;
    std::cout<<std::left<<"Acc : "<<Distance(x2, y)/Distance(y, zero)<<std::endl;
}

void Solution::InitMatrix(double down, double mid, double up)
{
    int n = A.size();
    A = std::vector<std::vector<double>>(n , std::vector<double>(n));
    A[0][0] = mid;
    A[0][1] = up;
    A[n - 1][n - 2] = down;
    for(int i = 1; i < n - 1; ++i)
    {
        A[i][i - 1] = down;
        A[i][i] = mid;
        A[i][i + 1] = up;
    }
    A[n - 1][n - 1] = mid;
}

void Solution::InitVector(double value)
{
    int n = b.size();
    b = std::vector<double>(n);
    for(int i = 0; i < n; ++i)
    {
        b[i] = value;
    }
}

double Solution::Distance(std::vector<double> x1, std::vector<double> x0)
{
    double distance = 0;

    for (int i=0;i<100;i++)
    {
        distance += (x1[i] - x0[i])*(x1[i] - x0[i]);
    }
    distance = sqrt(distance);

    return distance;
}

void Solution::AccuRes(double epsilon)
{
    for (int i=0;i<100;i++)
    {
        y[i] = 0.5 * (1 - exp(-0.01*i/epsilon)) /(1 - exp(-1/epsilon)) + 0.5*0.01*i;
    }
}

bool Solution::JudgeIteration(std::vector<double> x1, std::vector<double> x0, double e)
{
    int n = x1.size();
    double max_num = 0;
    for(int i = 0; i < n; ++i)
    {
        max_num = std::max(max_num, abs(x1[i] - x0[i]));
    }
    if(max_num <= e)
    {
        return false;
    }
    return true;
}