#include "utils.h"

string intToString(int value) {
    static char buf[64];
    sprintf(buf, "%d", value);
    return string(buf);
}

string doubleToString(double value) {
    static char buf[64];
    sprintf(buf, "%.02lf", value);
    return string(buf);
}
