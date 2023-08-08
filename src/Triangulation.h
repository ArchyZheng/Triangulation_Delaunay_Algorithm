//
// Created by 郑成琦 on 3/8/2023.
//

#ifndef LINEARINTERPOLATION_TRIANGULATION_H
#define LINEARINTERPOLATION_TRIANGULATION_H

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <set>
#include <list>

struct XYZ {
    std::string id; // the index of point according _samplePoints
    double x, y, z;
};

struct Triangle {
    std::set<std::string> id; // the index set of the points
    std::vector<XYZ> trianglePoints; // #trianglePoints = 3
    bool available = true; // if this triangle has been removed, set false.
};

enum StateOfPoint {
    outside, inside, onTheCircle
};

StateOfPoint isInsideTheTriangle(const XYZ &point, std::vector<XYZ> triangle);

bool checkTriangleContainSuperTriangleVertex(const Triangle &triangle);

void readFromBinaryFile(const std::string &fileStream, std::vector<double> *storedVector);

bool checkCross(XYZ *lineSegments1, XYZ *lineSegments2);

double getAngle(XYZ *vector1, XYZ *vector2);


/*
 * triangulation by Delaunay algorithm, ref: https://www.ics.uci.edu/~goodrich/teach/geom/notes/DT.pdf page:15
 */
class Triangulation {
public:
    Triangulation(std::vector<double> samplePointX, std::vector<double> samplePointY,
                  std::vector<double> samplePointZ);

    std::vector<Triangle> splitOneTriangleIntoThreeStateInside(const XYZ &point, std::vector<Triangle *> triangles);

    bool isTriangleCandidateContain(const Triangle& triangle);

    void splitOneTriangleIntoThreeStateOnTheCircle(XYZ point, Triangle triangle);

private:
    std::vector<XYZ> _samplePoints;
    std::list<Triangle> _triangleCandidate;
    Triangle _superTriangle;
};


#endif //LINEARINTERPOLATION_TRIANGULATION_H
