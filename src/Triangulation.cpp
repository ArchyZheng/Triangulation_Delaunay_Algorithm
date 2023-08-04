//
// Created by 郑成琦 on 3/8/2023.
//

#include "Triangulation.h"
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

/*
 *  estimate the sample point is inside, outside or on the circumcircle.
 *  getting the circumcircle center, and compare the distance between sample point and center.
 *  @parameter point: sample point
 *  @parameter triangle:
 */
StateOfPoint isInsideTheTriangle(XYZ point, std::vector<XYZ> triangle) {
    assert(triangle.size() == 3);
    const double EPSILON = 0.000001;

    double middleAB[] = {(triangle[0].x - triangle[1].x) / 2, (triangle[0].y - triangle[1].y) / 2};
    double middleAC[] = {(triangle[0].x - triangle[2].x) / 2, (triangle[0].y - triangle[2].y) / 2};
    // get the center of circumcircle.
    double slopeAB = (triangle[0].y - triangle[1].y) / (triangle[0].x - triangle[1].x + EPSILON);
    double slopeAC = (triangle[0].y - triangle[2].y) / (triangle[0].x - triangle[2].x + EPSILON);
    double slopeABMidperpendicular = -1 / (slopeAB + EPSILON);
    double slopeACMidperpendicular = -1 / (slopeAC + EPSILON);

    // the way of getting center location. ref: https://github.com/obviousjim/ofxDelaunay/blob/master/libs/Delaunay/src/Delaunay.cpp
    double centerX = -1 * (middleAB[0] * slopeABMidperpendicular - middleAC[0] * slopeACMidperpendicular + middleAC[1] -
                      middleAB[1]) / (slopeABMidperpendicular - slopeACMidperpendicular);
    double centerY = -1 * (slopeABMidperpendicular * (centerX - middleAB[0]) + middleAB[1]);

    double centerCircle[] = {centerX, centerY};
    auto distance = [](double *centerCircle, XYZ point) -> double {
        double distanceSquare = pow(centerCircle[0] - point.x, 2) + pow(centerCircle[1] - point.y, 2);
        return sqrt(distanceSquare);
    };

    // output part
    double redius = distance(centerCircle, triangle[0]);
    double distanceBetweenCenterAndSample = distance(centerCircle, point);
    if (redius > distanceBetweenCenterAndSample)
        return StateOfPoint::inside;
    if (redius < distanceBetweenCenterAndSample)
        return StateOfPoint::outside;
    return StateOfPoint::onTheCircle;
}
