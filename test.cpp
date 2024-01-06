#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include "utils.h"
#include <fstream>
#include <string>
#include <numeric>

TEST(vectorTest , Testemplace)
{
    std::vector<HSV*> origin_HSV(10);
    std::vector<int>  indice ;

    indice.push_back(3);  indice.push_back(5); indice.push_back(7);  
}

TEST(fileTest , Testwrite)
{
    std::string filename = "./models/bunny/reconstruction/bun_zipper_res3.ply";
    std::ifstream plyFile(filename);
    char delimiter = '/';

    std::vector<std::string> tokens = split(filename, delimiter);
    auto& ply_name =  *(tokens.end() - 1);
    char delimiter2 = '.';
    std::vector<std::string> ply_tokens = split(ply_name, delimiter2);
    auto ply_name2 =  *(ply_tokens.begin()) + "_seg.ply";
    std::string refine_filename =   ply_name2;
    std::cout << refine_filename << std::endl;


    std::string line;
    std::stringstream header;
    while (std::getline(plyFile, line)) {
        header << line << std::endl;
        if (line == "end_header") {
            break;
        }
    }

    std::ofstream refine_plyFile(refine_filename);

    refine_plyFile << header.str();
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}