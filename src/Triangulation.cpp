//
// Created by 郑成琦 on 3/8/2023.
//

#include "Triangulation.h"
#include "cmath"
#include <set>

/*
 * In construction function, it will generate a super triangle which involve every samplePoint.
 * The super triangle had been set by the previous knowledge which the location
 * of every sample point is inside [x in (-10, 10), y in (-10, 10)]
 *
 * @parameters:samplePointX the original points drawn triangle.
 * @parameters:samplePointY
 * @parameters:samplePointZ
 */
Triangulation::Triangulation(std::vector<double> samplePointX, std::vector<double> samplePointY,
                             std::vector<double> samplePointZ) {
    for (int i = 0; i < samplePointX.size(); i++) {
        std::string id = std::to_string(i);
        XYZ samplePoint = {id, samplePointX[i], samplePointY[i], samplePointZ[i]};
        this->_samplePoints.push_back(samplePoint);
    }

    // find a triangle which containing all points, 'super triangle'
    // we need three points, 1. (-30, -10) 2. (30, -10) 3. (0, 30)
    // adjust some location for getting each sample point inside the super triangle
    XYZ point1 = {"superTrianglePoint1", -30 - 30, -10 - 10, 0};
    XYZ point2 = {"superTrianglePoint2", 30 + 30, -10 - 10, 0};
    XYZ point3 = {"superTrianglePoint3", 0, 30 + 30, 0};

    std::vector<XYZ> superTrianglePoints = {point1, point2, point3};
    std::set < std::string > superTriangleId = {point1.id, point2.id, point3.id};
    Triangle superTriangle = {superTriangleId, superTrianglePoints};
    this->_triangleCandidate.push_back(superTriangle);
    this->_superTriangle = superTriangle;
    // repeat: select one point looking for triangle whose circumcircle contain this point
    // break those triangles and connect its edges to our point to make a new set of triangles.
    for (const auto &samplePoint: this->_samplePoints) {
        std::vector<Triangle> pointOnTheCircle; // the map will store triangles
        std::vector<Triangle> pointInsideTheCircle;
        for (const auto &triangle: this->_triangleCandidate) {
            StateOfPoint pointState = isInsideTheTriangle(samplePoint, triangle.trianglePoints);
            if (pointState == StateOfPoint::onTheCircle) {
                pointOnTheCircle.push_back(triangle);
            } else if (pointState == StateOfPoint::inside) {
                pointInsideTheCircle.push_back(triangle);
            }
        }
        for (auto triangle: pointOnTheCircle) {

        }
        for (const auto &triangle: pointInsideTheCircle) {
            splitOneTriangleIntoThreeStateInside(samplePoint, triangle);
        }
    }
    // delete all super-triangle vertices
}

/*
 * for the situation that the point is insider the triangle
 * splitting the triangle into three triangle by connect the point as the vertex of each triangle.
 * @parameter point: sample point
 * @parameter triangle: the point is inside the triangle
 * @return triangleVector: the output vector which contain three triangles
 */
std::vector<Triangle> Triangulation::splitOneTriangleIntoThreeStateInside(const XYZ &point, Triangle triangle) {
    std::vector<Triangle> triangleVectorOutput;
    int arrangement[3][2] = {{0, 1},
                             {0, 2},
                             {1, 2}};
    for (auto subArrangement: arrangement) {
        std::set < std::string > idSet = {point.id, triangle.trianglePoints[subArrangement[0]].id,
                                          triangle.trianglePoints[subArrangement[1]].id};
        std::vector<XYZ> pointVector = {point, triangle.trianglePoints[subArrangement[0]],
                                        triangle.trianglePoints[subArrangement[1]]};
        Triangle newTriangle = {idSet, pointVector};
        triangleVectorOutput.push_back(newTriangle);
    }
    return triangleVectorOutput;
}

/*
 *  estimate the sample point is inside, outside or on the circumcircle.
 *  getting the circumcircle center, and compare the distance between sample point and center.
 *  @parameter point: sample point
 *  @parameter triangle:
 */
StateOfPoint isInsideTheTriangle(const XYZ &point, std::vector<XYZ> triangle) {
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
    auto distance = [](double *centerCircle, const XYZ &point) -> double {
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

/*
 * read data from binary file, and restore it into vector outside.
 * read the data as double!
 *
 * @parameter filePath
 * @parameter storeVector
 */
void readFromBinaryFile(const std::string &filePath, std::vector<double> *storedVector) {
    std::ifstream fileStream(filePath, std::ios::binary);
    if (!fileStream.is_open()) {
        return;
    }
    fileStream.seekg(0, std::ios::end);
    auto fileSize = fileStream.tellg();
    fileStream.seekg(0, std::ios::beg);

    storedVector->resize(fileSize / sizeof(double));
    fileStream.read(reinterpret_cast<char *>(storedVector->data()), fileSize);
}
