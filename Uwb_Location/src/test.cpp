#include "stdio.h"

struct vec3d {
	double	x = 10000.0;
	double	y = 10000.0;
	double	z = 10000.0;
};

vec3d tmp;
printf("the primal point at(%f, %f, %f)\n", &tmp.x, &tmp.y, &tmp.z);