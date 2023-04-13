#include "func.h"



Solution::Solution(int n, std::vector<std::vector<double>> a1, std::string fn)
{
    if(n <= 0)    exit(1);
    A1.resize(n, std::vector<double>(n));
    A1 = a1;
    X.resize(n);
    Y.resize(n);
    b.resize(n);
    filename = fn;
}

void Solution::DoolittleDecom(int n)
{

    std::vector<std::vector<double>> A;
    
    A = A1;

    for(int i = 0; i < n; ++i)
    {
        for(int j = i; j < n; ++j)
        {
            U[i][j] = A[i][j];
            for(int k = 0; k < i; ++k)
            {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        for(int j = i; j < n; ++j)
        {
            L[j][i] = A[j][i];
            for(int k = 0; k < i; ++k)
            {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
    

}

void Solution::GetYofLY(int n, std::vector<std::vector<double>>& L1, std::vector<double>& b1, std::vector<double>& Y1)
{
    for(int i = 0; i < n; ++i)
    {
        Y1[i] = b1[i];
        for(int j = 0; j < i; ++j)
        {
            Y1[i] -= Y1[j] * L1[i][j];
        }
        Y1[i] = Y1[i] / L1[i][i];
    }
}

void Solution::GetXofUX(int n, std::vector<std::vector<double>>& U1, std::vector<double>& Y1, std::vector<double>& X1)
{
    for(int i = n - 1; i >= 0; --i)
    {
        X1[i] = Y1[i];
        for(int j = n - 1; j > i; --j)
        {
            X1[i] -= X1[j] * U1[i][j];
        }
        X1[i] /= U1[i][i];
    }
}

void Solution::MinEigIPM(int n, double epsilon)
{
    L.clear();
    L.resize(n, std::vector<double>(n));
    U.clear();
    U.resize(n, std::vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        X[i] = 1.0;
        L[i][i] = 1.0;
    }

    double lambda1 = 1.0;
    double lambda2 = 1.0;

    DoolittleDecom(n);
    int k = 0;

    std::ofstream outFile;
    outFile.open(filename, std::ios::out | std::ios::trunc);

    
    outFile<<"| k | Lambda |";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"$X_{"<<i + 1<<"}$ | ";
    }
    for(int i = 0; i < n; ++i)
    {
        outFile<<"$Y_{"<<i + 1<<"}$ | ";
    }
    outFile<<"\n";
    outFile<<"| ----------- | ----------- | ";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"----------- | ";
    }
    for(int i = 0; i < n; ++i)
    {
        outFile<<"----------- | ";
    }
    outFile<<"\n";
    outFile<<"0 | ";
    do
    {
        for (int i = 0; i < n; ++i)
        {
            Y[i] = X[i] / lambda2;
            outFile<<Y[i]<<" | ";
        }
        outFile<<"\n";
        lambda1 = lambda2;

        GetYofLY(n, L, Y, b);
        GetXofUX(n, U, b, X);

        lambda2 = X[0];
        for (int i = 1; i < n; ++i)
        {
            if (abs(X[i]) > abs(lambda2))
            {
                lambda2 = X[i];
            }
        }
        ++k;
        outFile<<k<<" | "<<lambda2<<" | ";
        for (int i=0;i<n;i++)   
        {
            outFile<<X[i]<<" | ";
        }
    } while (abs(lambda2 - lambda1) > epsilon);

    outFile<<"\n";
    outFile<<"\n";
    outFile<<"\n";
    outFile<<"Minimum eigenvalue is : "<<1/lambda2<<"\n\n\nEigenvector is :\n("<<Y[0];
    for (int i = 1;i < n; ++i)
    {
        outFile<<","<<Y[i];
    }
    outFile<<")\n";
    std::cout<<"Minimum eigenvalue is : "<<1/lambda2<<"\nEigenvector is :\n("<<Y[0];
    for (int i = 1; i < n; ++i)
    {
        std::cout<<", "<<Y[i];
    }
    std::cout<<")\n";
    outFile.close(); 

}