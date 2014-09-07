#ifndef GEOMUTILS_H
#define GEOMUTILS_H

#include <vector>
#include <math.h>
using namespace std;

#define eps 1e-7
const double Pi = 2.0 * acos(0.0);

struct pt {
    double x, y;
    pt() : x(0), y(0) {}
    pt(double X, double Y) : x(X), y(Y) {}
    double len();
    double len2();
    pt normalize();
};

bool eq(double a, double b);
bool lo(double a, double b);
double deg(double ang);
double rad(double ang);
pt operator +(pt a, pt b);
pt operator -(pt a, pt b);
pt operator *(pt a, double k);
double operator *(pt a, pt b);
double operator ^(pt a, pt b);
bool operator <(pt a, pt b);
bool operator ==(pt a, pt b);

struct circle {
    pt c;
    double r;
    circle(double x, double y, double R) { c = pt(x, y); r = R; }
    circle(pt p, double R) { c = p; r = R; }
};

struct line {
    double a, b, c;
    pt p1, p2;
    line() {}
    line(pt p1, pt p2) {
        this->p1 = p1;
        this->p2 = p2;
        a = p1.y-p2.y;
        b = p2.x-p1.x;
        c = -(a*p1.x+b*p1.y);
    }
};

// Line-to-circle intersection
vector<pt> line_circle_intersect(line l, circle c);

// Line-to-line intersection
bool line_line_intersect(line l1, line l2, pt& res);

// Get polar angle in the range [0, 2Pi) or -1 for (0, 0)
double polarAngle(pt a);

// Get rotation angle (ccw direction)
double rotationAngle(pt a, pt b);

// Rotate point t around c by angle alpha (ccw direction)
pt rotate(pt p, double alpha, pt c);

#endif // GEOMUTILS_H
