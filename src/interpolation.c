#include "interpolation.h"

float linear_int(float y0, float yh, float yth, float h){
    /*
    Uses linear interpolation of form y(t) = (yh-y0)/h * t + y0 to determine the intersection of y(t) with yth.
    */
   return h*(y0 - yth) / (y0 - yh);
}

float quadratic_int(float y0, float yh, float yth, float h);
float cubic_int(float y0, float yh, float yth, float h);