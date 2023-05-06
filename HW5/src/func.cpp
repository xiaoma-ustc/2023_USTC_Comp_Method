#include "func.h"

Solution::Solution()
{
}


void Solution::Thomas(int n, std::vector<std::vector<double>>& Mn3, std::vector<double>& b, std::vector<double>& x)
{
    double u = 0;
    std::vector<double> v(n);
    std::vector<double> y(n);
    y[0] = b[0] / Mn3[0][1];

    for(int i = 0; i < n; ++i)
    {
        u = Mn3[i][1] - Mn3[i][0] * v[i - 1];
        v[i] = Mn3[i][2] / u;
        y[i] = (b[i] - Mn3[i][0] * y[i - 1]) / u;
    }
    x[n - 1] = y[n - 1];
    for(int i = n - 2; i >= 0; --i)
    {
        x[i] = y[i] - v[i] * x[i + 1];
    }

}

void Solution::Spline(std::string FileName, std::vector<double>& Ms, int condition, int flag)
{
    std::vector<int> x;
    std::vector<double> y;
    std::ifstream inFile;
    inFile.open(FileName);
    for(int i = 0; i < 21; ++i)
    {
        double t1, t2;
        inFile>>t1>>t2;
        x.push_back(t1);
        y.push_back(t2);
    }
    inFile.close();
    int n = x.size();
    if(flag == 1)
    {
        x[9] = 0;
        y[9] = 10;
    }

    std::vector<std::vector<double>> A(n - 2, std::vector<double>(3));
    std::vector<double> M(n);
    std::vector<double> xM(n - 2);
    std::vector<double> b(n - 2);
    std::vector<double> mu(n - 1);
    std::vector<double> lambda(n - 1);
    std::vector<double> h(n - 1);
    std::vector<double> d(n - 1);

    h[0] = x[1] - x[0];
    for(int i = 1; i < (n - 1); ++i)
    {
        h[i] = x[i + 1] - x[i];
        lambda[i] = h[i] / (h[i] + h[i - 1]);
        mu[i] = 1 - lambda[i];
        d[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) / (h[i] + h[i - 1]);
    }

    if(condition == 0)
    {
        for(int i = 1; i < (n - 1); ++i)
        {
            A[i - 1][0] = mu[i];
            A[i - 1][1] = 2;
            A[i - 1][2] = lambda[i];
            b[i - 1] = d[i];
        }

    }
    Thomas(n - 2, A, b, xM);

    M[0] = 0;
    for(int i = 1; i < (n - 1); ++i)
    {
        M[i] = xM[i - 1];
    }
    M[n - 1] = 0;

    std::ofstream outFile;
    outFile.open("result" + std::to_string(flag) + ".csv", std::ios::out | std::ios::trunc);

    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    double a4 = 0;

    std::cout<<"\n";
    for(int i = 0; i < (n - 1); ++i)
    {
        a1 = (M[i+1] - M[i]) / ( 6 * h[i]);
        a2 = (x[i+1] * M[i] - x[i] * M[i+1]) / (2 * h[i]);
        a3 = (3 * (x[i] * x[i] * M[i + 1] - x[i + 1] * x[i + 1] * M[i])+ 6 * (y[i+1] - y[i]) - h[i] * h[i] * (M[i + 1] - M[i])) / (6 * h[i]);
        a4 = (x[i + 1] * x[i + 1] * x[i + 1] * M[i] - x[i] * x[i] * x[i] * M[i + 1] + 6 * (x[i + 1] * y[i] - x[i] * y[i + 1]) - h[i] * h[i] * (x[i + 1] * M[i] - x[i] * M[i + 1])) / (6 * h[i]);
        std::cout<<"["<<x[i]<<","<<x[i+1]<<"]  : ";
        std::cout<<"("<<a1<<")*x^3 + ";
        std::cout<<"("<<a2<<")*x^2 + ";
        std::cout<<"("<<a3<<")*x + ";
        std::cout<<"("<<a4<<")\n";

        outFile<<a1<<","<<a2<<","<<a3<<","<<a4<<"\n";
    }
    std::cout<<std::endl;

    outFile.close();
}