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

TEST(TriangulationSuit, splitSateInsideTriangleOnePoint) {
    XYZ point = {"samplePoint", 2, 1.5, 0};
    std::vector<double> samplePointX = {point.x};
    std::vector<double> samplePointY = {point.y};
    std::vector<double> samplePointZ = {point.z};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "hello world!" << std::endl;
}

TEST(TriangulationSuit, splitSateInsideTriangleThreePoints) {
    XYZ point1 = {"samplePoint", 1, 0, 0};
    XYZ point2 = {"samplePoint", 0, 0, 0};
    XYZ point3 = {"samplePoint", 0, 1, 0};
    std::vector<double> samplePointX = {point1.x, point2.x, point3.x};
    std::vector<double> samplePointY = {point1.y, point2.y, point3.y};
    std::vector<double> samplePointZ = {point1.z, point2.z, point3.z};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "hello world!" << std::endl;
}

TEST(TriangulationSuit, splitSateInsideTriangleFourPoints) {
    XYZ point1 = {"samplePoint", 1, 0, 0};
    XYZ point2 = {"samplePoint", 0, 0, 0};
    XYZ point3 = {"samplePoint", 0, 1, 0};
    XYZ point4 = {"samplePoint", 1, 1, 0};
    std::vector<double> samplePointX = {point1.x, point2.x, point3.x, point4.x};
    std::vector<double> samplePointY = {point1.y, point2.y, point3.y, point4.y};
    std::vector<double> samplePointZ = {point1.z, point2.z, point3.z, point4.z};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "hello world!" << std::endl;
}


TEST(TriangulationSuit, splitSateInsideTriangleFivePoints) {
    XYZ point1 = {"samplePoint", 1, 0, 0};
    XYZ point2 = {"samplePoint", 0, 0, 0};
    XYZ point3 = {"samplePoint", 0, 1, 0};
    XYZ point4 = {"samplePoint", 1, 1, 0};
    XYZ point5 = {"samplePoint", .5, -0.5, 0};
    std::vector<double> samplePointX = {point1.x, point2.x, point3.x, point4.x, point5.x};
    std::vector<double> samplePointY = {point1.y, point2.y, point3.y, point4.y, point5.y};
    std::vector<double> samplePointZ = {point1.z, point2.z, point3.z, point4.z, point5.z};
    Triangulation triangulation(samplePointX, samplePointY, samplePointZ);
    std::cout << "hello world!" << std::endl;
}

TEST(TriangulationSuit, checkCross) {
    XYZ point1 = {"samplePoint", 1, 0, 0};
    XYZ point2 = {"samplePoint", 0, 0, 0};
    XYZ point3 = {"samplePoint", .5, -1, 0};
    XYZ point4 = {"samplePoint", .5, 1, 0};
    XYZ lineOne[] = {point1, point2};
    XYZ lineTwo[] = {point3, point4};
    bool result = checkCross(lineOne, lineTwo);
    ASSERT_EQ(result, true);
    point1 = {"samplePoint", 1, 0, 0};
    point2 = {"samplePoint", 0, 0, 0};
    point3 = {"samplePoint", .5, -1, 0};
    point4 = {"samplePoint", 1, -1, 0};
    XYZ lineThree[] = {point1, point2};
    XYZ lineFour[] = {point3, point4};
    result = checkCross(lineThree, lineFour);
    ASSERT_EQ(result, false);

}
