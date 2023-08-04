//
// Created by 郑成琦 on 3/8/2023.
//
#include <iostream>
#include "Triangulation.h"
#include "gtest/gtest.h"

/*
 * test 1: super triangle
 */
TEST(TriangulationSuit, Construction) {
    std::vector<double> samplePointX = {10, 20 ,30, -40};
    std::vector<double> samplePointY = {-10, 20 ,-30, 40};
    std::vector<double> samplePointZ = {10, 20 ,30, 40};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "Test is running correct!" << std::endl;
}