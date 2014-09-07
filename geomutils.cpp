#include "geomutils.h"
#include <algorithm>
#include <iostream>

bool eq(double a, double b) { return fabs(a-b) < eps; }
bool lo(double a, double b) { return a+eps < b; }
double deg(double ang) { return ang*180/Pi; }
double rad(double ang) { return ang*Pi/180; }
pt operator +(pt a, pt b) { return pt(a.x+b.x, a.y+b.y); }
pt operator -(pt a, pt b) { return pt(a.x-b.x, a.y-b.y); }
pt operator *(pt a, double k) { return pt(k*a.x, k*a.y); }
double operator *(pt a, pt b) { return a.x*b.x + a.y*b.y; }
double operator ^(pt a, pt b) { return a.x*b.y - a.y*b.x; }
bool operator <(pt a, pt b) { return lo(a.x, b.x) || (eq(a.x, b.x) && lo(a.y, b.y)); }
bool operator ==(pt a, pt b) { return eq(a.x, b.x) && eq(a.y, b.y); }

double pt::len() { return sqrt(len2()); }
double pt::len2() { return x*x+y*y; }
pt pt::normalize() {
    double l = len();
    return l < eps ? pt(0, 0) : pt(x/l, y/l);
}

// Assume circle is placed in point (0, 0)
vector<pt> line_circle_intersect_simple(double r, double a, double b, double c) {
    vector<pt> ret;
    double t = a*a+b*b;
    double x0 = -a*c/t;
    double y0 = -b*c/t;
    if (c*c > r*r*t + eps)
        return ret; // No points
    else if (fabs(c*c-r*r*t) < eps)
        ret.push_back(pt(x0, y0)); // 1 point
    else { // 2 points
       double d = r*r - c*c/t;
       double mult = sqrt(d/t);
       ret.push_back(pt(x0+b*mult, y0-a*mult));
       ret.push_back(pt(x0-b*mult, y0+a*mult));
    }
    return ret;
}

vector<pt> line_circle_intersect(line l, circle c) {
    /*line L(pt(-5, -5), pt(5, 5));
    vector<pt> res = line_circle_intersect_simple(5, L.a, L.b, L.c);
    for (int i = 0; i < res.size(); i++) {
        cout << "{" << res[i].x << ";" << res[i].y << "}" << endl;
    }
    cout << "===" << endl;*/
    line l2(l.p1-c.c, l.p2-c.c);
    vector<pt> ret = line_circle_intersect_simple(c.r, l2.a, l2.b, l2.c);
    for (int i = 0; i < (int) ret.size(); i++) {
        ret[i].x += c.c.x;
        ret[i].y += c.c.y;
    }
    sort(ret.begin(), ret.end());
    return ret;
}

bool line_line_intersect(line l1, line l2, pt& res) {
    double d = l1.a*l2.b-l1.b*l2.a;
    if (fabs(d) < eps)
        return false;
    res.x = -(l1.c*l2.b-l1.b*l2.c)/d;
    res.y = -(l1.a*l2.c-l1.c*l2.a)/d;
    return true;
}

double polarAngle(pt a) {
   if (eq(a.x, 0) && eq(a.y, 0)) return -1;
   double p = atan2(a.y, a.x);
   return p < 0 ? p+2*Pi : p;
}

double rotationAngle(pt a, pt b) {
    double angA = polarAngle(a);
    double angB = polarAngle(b);
    double ret = angB - angA;
    if (ret < 0) ret += 2*Pi;
    return ret;
}

pt rotate(pt p, double alpha, pt c) {
   double acos = cos(alpha), asin = sin(alpha);
   return pt((p.x-c.x)*acos - (p.y-c.y)*asin + c.x,
             (p.x-c.x)*asin + (p.y-c.y)*acos + c.y);
}
