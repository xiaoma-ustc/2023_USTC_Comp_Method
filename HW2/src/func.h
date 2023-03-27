#ifndef FUNC
#define FUNC

#include<algorithm>
#include<vector>
class Solution
{
private:
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> y;
    std::vector<double> zero;
public:
    Solution(int n);
    ~Solution(){};

    void GaussEliminationWithPivoting(double epsilon);
    void GaussSeidelIteration(double epsilon);
    bool JudgeIteration(std::vector<double> x1, std::vector<double> x0, double e);
    void InitMatrix(double down, double mid, double up);
    void InitVector( double value);
    double Distance(std::vector<double> x1, std::vector<double> x0);
    void AccuRes(double epsilon);
};

#endif