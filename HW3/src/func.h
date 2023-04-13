#ifndef FUNC
#define FUNC
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <string>
class Solution
{
private:
    std::vector<std::vector<double>> A1;
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;
    std::vector<double> Y;
    std::vector<double> X;
    std::vector<double> b;
    std::string filename;
public:
    Solution(int n, std::vector<std::vector<double>> a1, std::string fn);
    void DoolittleDecom(int n);
    void MinEigIPM(int n, double epsilon);
    void GetYofLY(int n, std::vector<std::vector<double>>& L1, std::vector<double>& b1, std::vector<double>& Y1);
    void GetXofUX(int n, std::vector<std::vector<double>>& U1, std::vector<double>& Y1, std::vector<double>& X1);
};



#endif