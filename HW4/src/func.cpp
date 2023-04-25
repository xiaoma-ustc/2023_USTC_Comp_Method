#include "func.h"

Solution::Solution(int a, int b, std::vector<std::vector<double>>& iris)
{
    std::default_random_engine e;
    e.seed(time(0));
    A.resize(a, std::vector<double> (b));
    U.resize(a, std::vector<double> (a));
    Sigma.resize(a, std::vector<double> (a));
    VT.resize(b, std::vector<double> (b));
    for(int i = 0; i < a; ++i)
    {
        for(int j = 0; j < b; ++j)
        {
            
            std::uniform_real_distribution<double> u(0,1.0);
            A[i][j] = u(e);
        }
    }

    AT.resize(b, std::vector<double> (a));
    for(int i = 0; i < b; ++i)
    {
        for(int j = 0; j < a; ++j)
        {
            AT[i][j] = A[j][i];
        }
    }
    M = iris;
}

std::vector<std::vector<double>> Solution::returnA()
{
    return A;
}



std::vector<std::vector<double>> Solution::MultiMatrix(std::vector<std::vector<double>> X, std::vector<std::vector<double>> Y)
{
    int n = X.size();
    int m = X[0].size();
    std::vector<std::vector<double>> temp(n, std::vector<double> (n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            for(int k = 0; k < m; ++k)
            {
                temp[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
    return temp;
}

double Solution::VecDis(std::vector<double> X, std::vector<double> Y)
{
    int n = X.size();
    double temp;
    for(int i = 0; i < n; ++i)
    {
        temp += (X[i] - Y[i]) * (X[i] - Y[i]);
    }
    temp = sqrt(temp);
    return temp;
}

void Solution::VecNorm(std::vector<double>& X)
{
    double temp = 0;
    int n = X.size();
    for(int i = 0; i < n; ++i)
    {
        temp += X[i] * X[i];
    }
    temp = sqrt(temp);
    for(int i = 0; i < n; ++i)
    {
        X[i] /= temp;
    }
}

void Solution::GaussElimPartialPivot(std::vector<std::vector<double>>& X, std::vector<double>& Y, std::vector<double>& Z, int DetIs0)
{
    int n = X.size();
    for(int i = 0; i < n; ++i)
    {
        int iMax = i;
        for(int k = i + 1; k < n; ++k)
        {
            if(abs(X[iMax][i]) < abs(X[k][i]))
            {
                iMax = k;
            }
        }
        for(int j = 1; j < n; ++j)
        {
            double term;
            term = X[i][j];
            X[i][j] = X[iMax][j];
            X[iMax][j] = term;
        }
        double term = Z[i];
        Z[i] = Z[iMax];
        Z[iMax] = term;

        for(int k = i + 1; k < n; ++k)
        {
            X[k][i] = X[k][i] / X[i][i];
            for(int j = i + 1; j < n; ++j)
            {
                X[k][j] -= X[k][i] * X[i][j];
            }
            Z[k] -= X[k][i] * Z[i];
        }
    }

    if(DetIs0 != -1 && DetIs0 != 0 && DetIs0 != 1)
    {
        exit(1);
    }

    if(DetIs0 == -1)
    {
        double LambdaMin = X[0][0];
        double LambdaMax = X[0][0];
        for(int k = 0; k < n; ++k)
        {
            if(abs(X[k][k]) > LambdaMax)
            {
                LambdaMax = X[k][k];
            }
            if(abs(X[k][k]) < LambdaMin)
            {
                LambdaMin = X[k][k];
            }
        }
        DetIs0 = (LambdaMin < 1 && abs(LambdaMin / LambdaMax) < 0.005) ? 1 : 0;
    }

    if(DetIs0 == 0)
    {
        for(int i = n - 1; i >= 0; --i)
        {
            for(int j = i + 1; j < n; ++j)
            {
                Z[i] -= Z[j] * X[i][j];
            }
            Z[i] /= X[i][i];
        }
    }
    else
    {
        if(Z[n - 1] != 0.0)
        {
            std::cout<<"***No Solution!***"<<std::endl;
            exit(1);
        }
        Z[n - 1] = 1.0;
        for(int i = n - 1; i >= 0; --i)
        {
            for(int j = i + 1; j < n; ++j)
            {
                Z[i] -= Z[j] * X[i][j];
            }
            Z[i] /= X[i][i];
        }
        VecNorm(Z);
    }
    for(int i = 0; i < n; ++i)
    {
        Y[i] = Z[i];
    }
}


void Solution::GaussSeidel(std::vector<std::vector<double>>& X, std::vector<double>& Y, std::vector<double>& Z, double epsilon)
{
    int n = X.size();
    std::vector<double> x1(n, 1);
    std::vector<double> x0(n, 10);
    std::vector<std::vector<double>> M(n, std::vector<double> (n));
    
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            M[i][j] = (i == j? 0 : (-X[i][j] / X[i][i]));
        }
    }

    for(int i = 0; i < n; ++i)
    {
        Z[i] /= X[i][i];
    }

    while(VecDis(x1, x0) > epsilon)
    {
        for(int i = 0; i < n; ++i)
        {
            x0[i] = x1[i];
        }
        for(int i = 0; i < n; ++i)
        {
            x1[i] = Z[i];
            for(int j = 0; j < n; ++j)
            {
                x1[i] += M[i][j] * x1[j];
            }
        }
    }

    for(int i = 0; i < n; ++i)
    {
        Y[i] = x1[i];
    }

}


void Solution::EigValJacobi(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVal, double epsilon, int flag)
{
    

    double sum_diag2 = 0;
    double MaxDiag = 0;
    int MaxDiag_p;
    int MaxDiag_q;
    double s, t, cos, sin;
    int n = Matrix.size();
    std::vector<std::vector<double>> temp = Matrix;

    
    do
    {
        sum_diag2 = 0.0;
        MaxDiag = 0.0;
        MaxDiag_p = 0;
        MaxDiag_q = 0;
        
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                if((i != j) && (fabs(MaxDiag) < fabs(temp[i][j])))
                {
                    MaxDiag = temp[i][j];
                    MaxDiag_p = i;
                    MaxDiag_q = j;
                }
            }
        }

        s = (temp[MaxDiag_q][MaxDiag_q] - temp[MaxDiag_p][MaxDiag_p]) / (2.0 * temp[MaxDiag_p][MaxDiag_q]);
        
        if(s == 0)
        {
            t = 1.0;
            cos = 1.0 / sqrt(2);
            sin = 1.0 / sqrt(2);
        }
        else
        {
            
            t = (s > 0 ? (-s + sqrt(1 + s * s)) : (-s - sqrt(1 + s * s)));
            
            cos = 1.0 / sqrt(1 + t * t);
            sin = t / sqrt(1 + t * t);
        }
        
        for(int i = 0; i < n; ++i)
        {
            if((i == MaxDiag_p) || (i == MaxDiag_q))
            {
                continue;
            }
            double B_ip = temp[MaxDiag_p][i] * cos - temp[MaxDiag_q][i] * sin;
            double B_iq = temp[MaxDiag_p][i] * sin + temp[MaxDiag_q][i] * cos;

            temp[MaxDiag_p][i] = B_ip;
            temp[i][MaxDiag_p] = B_ip;
            temp[MaxDiag_q][i] = B_iq;
            temp[i][MaxDiag_q] = B_iq;


        }
        double B_pp = temp[MaxDiag_p][MaxDiag_p] * cos * cos + temp[MaxDiag_q][MaxDiag_q] * sin * sin - temp[MaxDiag_p][MaxDiag_q] * 2.0 * sin * cos;
        double B_qq = temp[MaxDiag_p][MaxDiag_p] * sin * sin + temp[MaxDiag_q][MaxDiag_q] * cos * cos + temp[MaxDiag_p][MaxDiag_q] * 2.0 * sin * cos;
        double B_pq = temp[MaxDiag_p][MaxDiag_q] * (cos * cos - sin * sin) + (temp[MaxDiag_p][MaxDiag_p] - temp[MaxDiag_q][MaxDiag_q]) * sin * cos;
        
        temp[MaxDiag_p][MaxDiag_p] = B_pp;
        temp[MaxDiag_q][MaxDiag_q] = B_qq;
        temp[MaxDiag_q][MaxDiag_p] = B_pq;
        temp[MaxDiag_p][MaxDiag_q] = B_pq;
        
        
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                sum_diag2 += (i == j? 0 : temp[i][j] * temp[i][j]);
            }
        }
        if(flag == 1)
            sum.push_back(sum_diag2);

        
    }while(sum_diag2 > epsilon);


    
    //std::cout<<"EigVal :"<<std::endl;
    for(int i = 0; i < n; ++i)
    {
        EigVal[i] = temp[i][i];
        //std::cout<<EigVal[i]<<" ";
    }
    //std::cout<<std::endl; 



}

void Solution::EigVec(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVec, double EigVal)
{
    int n = Matrix.size();
    std::vector<double> b(n);
    std::vector<std::vector<double>> temp(n, std::vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            temp[i][j] = (i == j? Matrix[i][j] - EigVal : Matrix[i][j]);
        }
    }

    GaussElimPartialPivot(temp, EigVec, b, 1);

}

/* void Solution::EigVecIter(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVec, double EigVal, double epsilon)
{
    int n = Matrix.size();
    std::vector<double> b(n);
    std::vector<std::vector<double>> temp(n, std::vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            temp[i][j] = (i == j? Matrix[i][j] - EigVal : Matrix[i][j]);
        }
    }

    GaussElimPartialPivot(temp, EigVec, b, epsilon);
    
} */

void Solution::SVD(double epsilon)
{
    int n = A.size();
    int m = A[0].size();
    std::vector<double> AAT_EigVal(n);
    std::vector<double> ATA_EigVal(m);
    AAT = MultiMatrix(A, AT);
    ATA = MultiMatrix(AT, A);

    

    EigValJacobi(AAT, AAT_EigVal, epsilon, 1);
    
    std::sort(AAT_EigVal.begin(), AAT_EigVal.end(), std::greater<double>());
    
    for(int i = 0; i < n; ++i)
    {
        std::vector<double> eigvec(n);
        EigVec(AAT, eigvec, AAT_EigVal[i]);
        for(int j = 0; j < n; ++j)
        {
            U[j][i] = eigvec[j];
        }
    }

    EigValJacobi(ATA, ATA_EigVal, epsilon,0);
    std::sort(ATA_EigVal.begin(), ATA_EigVal.end(), std::greater<double>());
    
    for(int i = 0; i < m; ++i)
    {
        std::vector<double> eigvec(n);
        EigVec(ATA, eigvec, ATA_EigVal[i]);
        for(int j = 0; j < m; ++j)
        {
            VT[i][j] = eigvec[j];
        }
    }

    //Sigma.resize(n, std::vector<double> (m));
    //double min_n_m = m < n? m : n;
    for(int i = 0; i < m; ++i)
    {
        Sigma[i][i] = sqrt(AAT_EigVal[i]);
        
    }


}

void Solution::PCA(double epsilon, int FinalDim, std::vector<std::vector<double>>& FinalM, int Buffer)
{
    int n = M.size();
    int m = M[0].size();
    std::vector<std::vector<double>> temp = M;
    

    std::vector<double> meanvec(n);
    for(int j = 0; j < m; ++j)
    {
        for(int i = 0; i < n; ++i)
        {
            meanvec[i] += temp[i][j];
        }
    }
    for(int i = 0; i < n; ++i)
    {
        meanvec[i] /= m;
    }
    for(int j = 0; j < m; ++j)
    {
        for(int i = 0; i < n; ++i)
        {
            temp[i][j] -= meanvec[i];
        }
    }

    std::vector<std::vector<double>> MT(m, std::vector<double>(n));
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            MT[i][j] = temp[j][i];
        }
    }
    std::vector<std::vector<double>> MMT;
    MMT = MultiMatrix(temp, MT);
    std::ofstream outFile;
    outFile.open("PCA.md", std::ios::out | std::ios::trunc);
    outFile<<"1/m XXT \n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            outFile<<"| "<<MMT[i][j]/m;
        }
        outFile<<"|\n";
    }
    outFile<<"\n";

    std::vector<double> MMT_EigVal(n);
    
    
    EigValJacobi(MMT, MMT_EigVal, epsilon,0);
    sort(MMT_EigVal.begin(), MMT_EigVal.end(), std::greater<double>());
    
    for(int i = 0; i < FinalDim; ++i)
    {
        std::vector<double> base(n);
        EigVec(MMT, base, MMT_EigVal[i]);
        
        for(int j = 0; j < m; ++j)
        {
            FinalM[i][j] = 0;
            for(int k = 0; k < n; ++k)
            {
                FinalM[i][j] += temp[k][j] * base[k];
            }
        }
        //std::cout<<FinalM[0][0]<<" "<<FinalM[1][0];
    }
}

void Solution::PrintSVD()
{
    int n = A.size();
    int m = A[0].size();
    std::ofstream outFile;
    outFile.open("SVD.md", std::ios::out | std::ios::trunc);
    outFile<<"A \n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
        {
            outFile<<"| "<<A[i][j];
        }
        outFile<<"|\n";
    }
    outFile<<"\n";
    outFile<<"AAT \n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            outFile<<"| "<<AAT[i][j];
        }
        outFile<<"|\n";
    }
    outFile<<"\n";
    outFile<<"U \n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            outFile<<"| "<<U[i][j];
        }
        outFile<<"|\n";
    }
    outFile<<"\n";
    outFile<<"Sigma \n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            outFile<<"| "<<Sigma[i][j];
        }
        outFile<<"|\n";
    }
    outFile<<"\n";
    outFile<<"VT \n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < m; ++i)
    {
        outFile<<"| ------ ";
    }
    outFile<<"|";
    outFile<<"\n";
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < m; ++j)
        {
            outFile<<"| "<<VT[i][j];
        }
        outFile<<"|\n";
    }
    outFile<<"\n";
    outFile<<"Sum_N_Diag\n";
    outFile<<"\n";
    for(auto d : sum)
    {
        outFile<<d<<"\n";
        outFile<<"\n";
    }
    
    outFile.close();
}