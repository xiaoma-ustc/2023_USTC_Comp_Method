#include "func.h"
#include "func.cpp"
#include <sstream>
int main()
{
    double epsilon = 0.000001;
    std::vector<std::vector<double>> iris(4, std::vector<double> (150));
    std::vector<std::vector<double>> FinalM(2, std::vector<double> (150));
    std::vector<int> label(150);
    std::ifstream inFile;
    inFile.open("iris.txt", std::ios::in);
    for(int i = 0; i < 150; ++i)
    {
        std::string s;
        if (!getline(inFile, s)) break;

        std::istringstream ss(s);
        std::vector <std::string> record;

        while (ss)
        {
            std::string s;
            if (!std::getline(ss, s, ',')) break;
            record.push_back(s);
            
        }

        iris[0][i] = atof(record[0].c_str());
        iris[1][i] = atof(record[1].c_str());
        iris[2][i] = atof(record[2].c_str());
        iris[3][i] = atof(record[3].c_str());
        label[i] = atof(record[4].c_str());;
        //std::cout<<iris[0][i]<<","<<iris[1][i]<<","<<iris[2][i]<<","<<iris[3][i]<<std::endl;
    }
    inFile.close();


    Solution sol = Solution(4, 3, iris);
    sol.SVD(epsilon);
    sol.PrintSVD();
    sol.PCA(epsilon, 2, FinalM, 1);
    std::ofstream outFile;
    outFile.open("result.csv", std::ios::out | std::ios::trunc);
    for(int i = 0; i < 150; ++i)
    {
        outFile<<FinalM[0][i]<<","<<FinalM[1][i]<<","<<label[i]<<"\n";
    }
    outFile.close();
    return 0;
}