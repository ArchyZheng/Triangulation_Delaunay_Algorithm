//
// Created by 郑成琦 on 3/8/2023.
//

#include "Triangulation.h"
#include "cmath"

/*
 * In construction function, it will generate a super triangle which involve every samplePoint.
 * The super triangle had been set by the previous knowledge which the location
 * of every sample point is [(-10, 10), (-10, 10)]
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
    // we need three points, 1. (-30, -10) 2. (30, -10) 3. (0, 30)
    // adjust some location for getting each sample point inside the super triangle
    XYZ point1 = {-30 - 30, -10 - 10, 0};
    XYZ point2 = {30 + 30, -10 - 10, 0};
    XYZ point3 = {0, 30 + 30, 0};

    std::vector<XYZ> superTriangle = {point1, point2, point3};
    this->_triangleCandidate.insert(std::pair<std::string, std::vector<XYZ>>("SuperTriangle", superTriangle));
    // repeat: select one point look for triangle whose circumcircle contain this point
    // break those triangles and connect its edges to our point to make a new set of triangles.
    for (auto samplePoint: this->_samplePoints) {
       std::map<std::string, std::vector<XYZ>> pointOnTheCircle; // the map will store triangles on which the point
        std::map<std::string, std::vector<XYZ>> pointInsideTheCircle; // the map will store triangles inside which the point
        for (const auto &triangle: this->_triangleCandidate) {
            StateOfPoint pointState = isInsideTheTriangle(samplePoint, triangle.second);
            if (pointState == StateOfPoint::onTheCircle) {
                pointOnTheCircle.insert(triangle);
            } else if (pointState == StateOfPoint::inside) {
                pointInsideTheCircle.insert(triangle);
            }
        }
    }
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
/*
 * read data from binary file, and restore it into vector outside.
 * read the data as double!
 *
 * @parameter filePath
 * @parameter storeVector
 */
void readFromBinaryFile(const std::string& filePath, std::vector<double> *storedVector) {
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
