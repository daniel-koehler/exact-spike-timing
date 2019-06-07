#ifndef interpolation_h
#define interpolation_h

float linear_int(float y0, float yh, float yth, float h);
float quadratic_int(float y0, float yh, float yth, float h);
float cubic_int(float y0, float yh, float yth, float h);

#endif