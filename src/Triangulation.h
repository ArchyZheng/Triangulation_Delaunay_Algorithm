//
// Created by 郑成琦 on 3/8/2023.
//

#ifndef LINEARINTERPOLATION_TRIANGULATION_H
#define LINEARINTERPOLATION_TRIANGULATION_H

#include <vector>
#include <map>
#include <string>

struct XYZ {
    double x, y, z;
};

enum StateOfPoint {
    outside, inside, onTheCircle
};

StateOfPoint isInsideTheTriangle(XYZ point, std::vector<XYZ> triangle);

class Triangulation {
public:
    Triangulation(std::vector<double> samplePointX, std::vector<double> samplePointY, std::vector<double> samplePointZ);

private:
    std::vector<XYZ> _samplePoints;
    std::map<std::string, std::vector<XYZ>> _triangleCandidate;
    std::vector<XYZ> _superTriangle;
};


#endif //LINEARINTERPOLATION_TRIANGULATION_H
