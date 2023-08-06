//
// Created by 郑成琦 on 3/8/2023.
//
#include <iostream>
#include "Triangulation.h"
#include "gtest/gtest.h"
#include "fstream"

/*
 * test 1: super triangle
 */
TEST(TriangulationSuit, Construction) {
    std::vector<double> samplePointX = {1.3, 1.5, 2.9, -3.5};
    std::vector<double> samplePointY = {9.6, -8.7, 6.7, -2.8};
    std::vector<double> samplePointZ = {10, 20, 30, 40};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "Test is running correct!" << std::endl;
}

TEST(TriangulationSuit, isInsideTheTriangle) {
    std::vector<XYZ> triangle = {{"point1", 0, 0, 0},
                                 {"point2", 0, 3, 0},
                                 {"point3", 4, 0, 0}};
    XYZ point = {"samplePoint", 2, 1.5, 0};
    ASSERT_EQ(isInsideTheTriangle(point = point, triangle = triangle), StateOfPoint::inside);
    point = {"samplePoint", 0, 0, 0};
    ASSERT_EQ(isInsideTheTriangle(point, triangle), StateOfPoint::onTheCircle);
    point = {"samplePoint", 100, 100, 0};
    ASSERT_EQ(isInsideTheTriangle(point, triangle), StateOfPoint::outside);
}

TEST(TriangulationSuit, loadData) {
    std::vector<double> sampleXVector;
    std::vector<double> sampleYVector;
    std::vector<double> sampleZVector;
    std::string sampleXPath = "../data/x_src.bin";
    std::string sampleYPath = "../data/y_src.bin";
    std::string sampleZPath = "../data/z_src.bin";
    readFromBinaryFile(sampleXPath, &sampleXVector);
    readFromBinaryFile(sampleYPath, &sampleYVector);
    readFromBinaryFile(sampleZPath, &sampleZVector);
    Triangulation triangulation(sampleXVector, sampleYVector, sampleZVector);
}
