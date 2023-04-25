#ifndef FUNC
#define FUNC

#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <fstream>
#include <algorithm>

class Solution
{
private:
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> AT;
    std::vector<std::vector<double>> AAT;
    std::vector<std::vector<double>> ATA;
    std::vector<std::vector<double>> U;
    std::vector<std::vector<double>> Sigma;
    std::vector<std::vector<double>> VT;
    std::vector<std::vector<double>> M;
    std::vector<double> sum;
public:
    Solution(int a, int b, std::vector<std::vector<double>>& iris);
    std::vector<std::vector<double>> returnA();
    std::vector<std::vector<double>> MultiMatrix(std::vector<std::vector<double>> X, std::vector<std::vector<double>> Y);
    double VecDis(std::vector<double> X, std::vector<double> Y);
    void GaussElimPartialPivot(std::vector<std::vector<double>>& X, std::vector<double>& Y, std::vector<double>& Z, int DetIs0);
    void VecNorm(std::vector<double>& X);
    void GaussSeidel(std::vector<std::vector<double>>& X, std::vector<double>& Y, std::vector<double>& Z, double epsilon);
    void EigValJacobi(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVal, double epsilon, int flag);
    void EigVec(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVec, double EigVal);
    //void EigVecIter(std::vector<std::vector<double>>& Matrix, std::vector<double>& EigVec, double EigVal, double epsilon);
    void SVD(double epsilon);
    void PCA(double epsilon, int FinalDim, std::vector<std::vector<double>>& FinalM, int Buffer);
    void PrintSVD();


};



#endif