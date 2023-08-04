//
// Created by 郑成琦 on 3/8/2023.
//

#include "Triangulation.h"
#include "iostream"
#include "cmath"

/*
 * In construction function, it will generate a super triangle which involve every samplePoint.
 *
 * @parameters:samplePointX the original point as the vertex of drawn triangle.
 * @parameters:samplePointY
 * @parameters:samplePointZ
 */
Triangulation::Triangulation(std::vector<double> samplePointX, std::vector<double> samplePointY,
                             std::vector<double> samplePointZ) {
    for (int i = 0; i < samplePointX.size(); i++) {
        XYZ samplePoint = {samplePointX[i], samplePointY[i], samplePointZ[i]};
        this->_samplePoints.push_back(samplePoint);
    }

    // find a triangle which containing all points, 'super triangle'
    // we need three points,
    // 1. the left most point, this point comes from samplePoint 2. right most and top 3.right and button.
    // last two point is created by hands.
    auto leftMostPoint = std::min_element(samplePointX.begin(), samplePointX.end());
    auto leftMostPointIndex = std::distance(samplePointX.begin(), leftMostPoint);
    XYZ point1 = {samplePointX[leftMostPointIndex], samplePointY[leftMostPointIndex], samplePointZ[leftMostPointIndex]};

    auto rightMostPoint = std::max_element(samplePointX.begin(), samplePointX.end());
    auto rightMostPointIndex = std::distance(samplePointX.begin(), rightMostPoint);

    // add 100 as bias, this will make the line between point1 and point2 will not pass any sample point.
    double rightAdjust = samplePointX[rightMostPointIndex] + 100;

    XYZ point2 = {rightAdjust, INFINITY, 0};
    XYZ point3 = {rightAdjust, -INFINITY, 0};
    std::vector<XYZ> superTriangle = {point1, point2, point3};
    this->_triangleCandidate.insert(std::pair<std::string, std::vector<XYZ>>("SuperTriangle", superTriangle));

    // repeat: select one point look for triangle whose circumcircle contain this point
    // break those triangles and connect its edges to our point to make a new set of triangles.

    // delete all super-triangle vertices
}
