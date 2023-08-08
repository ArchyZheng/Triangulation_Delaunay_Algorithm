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
        std::vector<Triangle *> pointOnTheCircle; // the map will store triangles
        std::vector<Triangle *> pointInsideTheCircle;
        for (auto &triangle: this->_triangleCandidate) {
            if (!triangle.available)
                continue;
            StateOfPoint pointState = isInsideTheTriangle(samplePoint, triangle.trianglePoints);
            if (pointState == StateOfPoint::onTheCircle) {
                pointOnTheCircle.push_back(&triangle);
            } else if (pointState == StateOfPoint::inside) {
                pointInsideTheCircle.push_back(&triangle);
            }
        }
        std::vector<Triangle> addedTriangle = splitOneTriangleIntoThreeStateInside(samplePoint,
                                                                                   pointInsideTheCircle);
        for (const auto &triangle: addedTriangle) {
            this->_triangleCandidate.push_back(triangle);
        }
    }
    // delete all super-triangle vertices
    std::list<Triangle> output;
    for (const auto &triangle: this->_triangleCandidate) {
        if (!triangle.available) {
            continue;
        }
        if (checkTriangleContainSuperTriangleVertex(triangle)) {
            continue;
        }
        output.push_back(triangle);
    }
    this->_triangleCandidate = output;
}

/*
 * this function is for clearing up the triangle vector
 * @return ture if this triangle contain superTriangleVertex. else return false
 */
bool checkTriangleContainSuperTriangleVertex(const Triangle &triangle) {
    std::vector<std::string> superTrianglePoints = {"superTrianglePoint1", "superTrianglePoint2",
                                                    "superTrianglePoint3"};
    for (const auto &superTrianglePoint: superTrianglePoints) {
        if (triangle.id.count(superTrianglePoint)) {
            return true;
        }
    }

    return false;
}

/*
 * Bowyer Watson algorithm ref: https://www.youtube.com/watch?v=4ySSsESzw2Y&t=256s.
 * point will connect each point of triangle, to create new triangles.
 * @parameter point: sample point
 * @parameter trianglesContainThePoint: triangles contain the point
 * @return triangleVector: the output vector which contain three triangles
 */
std::vector<Triangle>
Triangulation::splitOneTriangleIntoThreeStateInside(const XYZ &point,
                                                    std::vector<Triangle *> trianglesContainThePoint) {
    std::vector<std::pair<XYZ, double>> angleBetweenPointAndTrianglePoints;

    XYZ standardVector[] = {{"base", 1, 0, 0},
                            {"base", 0, 0, 0}}; // for compute the angle of each segmentation.
    std::set < std::string > attendPoint;
    for (auto &triangle: trianglesContainThePoint) {
        // make each triangle in the vector as bad triangle
        triangle->available = false;

        // get each angle, and sort it
        for (const auto &trianglePoint: triangle->trianglePoints) {
            XYZ vector[] = {point, trianglePoint};
            double angle = getAngle(vector, standardVector);
            std::string temp = trianglePoint.id;
            if (attendPoint.count(temp)) {
                continue;
            }
            attendPoint.insert(temp);
            auto storedAngle = std::make_pair(trianglePoint, angle);
            angleBetweenPointAndTrianglePoints.push_back(storedAngle);
        }
        std::sort(angleBetweenPointAndTrianglePoints.begin(), angleBetweenPointAndTrianglePoints.end(),
                  [](const auto &a, const auto &b) {
                      return a.second > b.second;
                  });
    }
    // point connects each point in the triangles
    angleBetweenPointAndTrianglePoints.push_back(
            angleBetweenPointAndTrianglePoints[0]); // connect the final triangle of the circle
    std::vector<Triangle> outputTriangle;
    for (int i = 0; i < angleBetweenPointAndTrianglePoints.size() - 1; i++) {
        // reassign the triangle list
        std::set < std::string > triangleId = {point.id, angleBetweenPointAndTrianglePoints[i].first.id,
                                               angleBetweenPointAndTrianglePoints[i + 1].first.id};
        std::vector<XYZ> trianglePoints = {point, angleBetweenPointAndTrianglePoints[i].first,
                                           angleBetweenPointAndTrianglePoints[i + 1].first};
        Triangle newTriangle = {triangleId, trianglePoints, true};
        outputTriangle.push_back(newTriangle);
    }
    return outputTriangle;
}


bool Triangulation::isTriangleCandidateContain(const Triangle &triangle) {
    for (const auto &triangleFromCandidate: this->_triangleCandidate) {
        if (triangleFromCandidate.id == triangle.id) {
            return true;
        }
    }
    return false;
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

/*
 * this function will check whether two line segmentation cross.
 * ref https://www.cnblogs.com/tuyang1129/p/9390376.html
 * @return ture if cross, else return false
 */
bool checkCross(XYZ *lineSegments1, XYZ *lineSegments2) {
    double vector1[2] = {lineSegments1[0].x - lineSegments1[1].x, lineSegments1[0].y - lineSegments1[1].y};
    double vector2[2] = {lineSegments2[0].x - lineSegments2[1].x, lineSegments2[0].y - lineSegments2[1].y};

    double vector1_1[2] = {lineSegments2[0].x - lineSegments1[0].x, lineSegments2[0].y - lineSegments1[0].y};
    double vector1_2[2] = {lineSegments2[1].x - lineSegments1[0].x, lineSegments2[1].y - lineSegments1[0].y};

    double vector2_1[2] = {lineSegments1[0].x - lineSegments2[0].x, lineSegments1[0].y - lineSegments2[0].y};
    double vector2_2[2] = {lineSegments1[1].x - lineSegments2[0].x, lineSegments1[1].y - lineSegments2[0].y};

    bool isTwoSideForFirstLine = (vector1[0] * vector1_1[1] - vector1[1] * vector1_1[0]) *
                                 (vector1[0] * vector1_2[1] - vector1[1] * vector1_2[0]) < 0;
    bool isTwoSideForSecondLine = (vector2[0] * vector2_1[1] - vector2[1] * vector2_1[0]) *
                                  (vector2[0] * vector2_2[1] - vector2[1] * vector2_2[0]) < 0;
    return isTwoSideForFirstLine && isTwoSideForSecondLine;
}

/*
 * as the function name.
 * @parameter vector1
 * @parameter vector2
 */
double getAngle(XYZ *vector1, XYZ *vector2) {
    double lineSegmentation1[] = {vector1[0].x - vector1[1].x, vector1[0].y - vector1[1].y};
    double lineSegmentation2[] = {vector2[0].x - vector2[1].x, vector2[0].y - vector2[1].y};

    double innerProduct = lineSegmentation1[0] * lineSegmentation2[0] + lineSegmentation1[1] * lineSegmentation2[1];

    double length1 = std::sqrt(pow(lineSegmentation1[0], 2.0) + pow(lineSegmentation1[1], 2.0));
    double length2 = std::sqrt(pow(lineSegmentation2[0], 2.0) + pow(lineSegmentation2[1], 2.0));

//    double direction = lineSegmentation1[0] * lineSegmentation2[1] - lineSegmentation1[1] * lineSegmentation2[0];
//    float rho = 1;
//    if (direction < 0)
//        rho = -1;

    return acos(innerProduct / (length1 * length2));
}
