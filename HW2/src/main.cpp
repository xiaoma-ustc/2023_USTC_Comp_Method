#include<iostream>
#include"func.h"
#include"func.cpp"
int main()
{//1.0,0.1,0.01,
    Solution sol = Solution(100);
    std::vector<double> epsilons{1.0,0.1,0.01,0.0001};
    double h = 0.01;
    for(auto epsilon : epsilons)
    {
        std::cout<<"epsilon :"<<epsilon<<std::endl;
        sol.AccuRes(epsilon);
        sol.InitMatrix(epsilon, -(2 * epsilon + h), epsilon + h);
        sol.GaussEliminationWithPivoting(epsilon);
        std::cout<<"epsilon :"<<epsilon<<std::endl;
        sol.InitMatrix(epsilon, -(2 * epsilon + h), epsilon + h);
        sol.GaussSeidelIteration(epsilon);

    }
    return 0;
}