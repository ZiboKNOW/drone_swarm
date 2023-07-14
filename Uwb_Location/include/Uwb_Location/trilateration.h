// -------------------------------------------------------------------------------------------------------------------
//
//  File: trilateration.h
//
//  Copyright 2016 (c) Decawave Ltd, Dublin, Ireland.
//
//  All rights reserved.
//
//
// -------------------------------------------------------------------------------------------------------------------

#ifndef __TRILATERATION_H__
#define __TRILATERATION_H__

#include "stdio.h"

#define		TRIL_3SPHERES							3
#define		TRIL_4SPHERES							4
#define   Kp                        2.0f
#define   Ki                        0.2f
#define   halfT                     0.05f
#define   squa( Sq )                (((float)Sq)*((float)Sq))
#define   MAX_AHCHOR_NUMBER         8
//四元数姿态解算
//float NormAccz;
//const float M_PI = 3.1415926535;
const float RtA = 57.2957795f;
const float Gyro_G = 0.03051756f;	//陀螺仪int16角速度除以分辨率得到度 量程1000
const float Gyro_Gr = 0.0005326f; //陀螺仪int16角速度转度再转弧度 量程1000
//static Quaternion NumQ = {1, 0, 0, 0};  // 四元素

float Q_rsqrt(float number);
typedef struct vec3d	vec3d;

struct vec3d {
	double	x = 10000.0;
	double	y = 10000.0;
	double	z = 10000.0;
};

//typedef char uint8_t;
typedef struct
{
    __int16_t accX;
    __int16_t accY;
    __int16_t accZ;
    __int16_t gyroX;
    __int16_t gyroY;
    __int16_t gyroZ;
}ImuData_t;

typedef struct
{
	float roll;
	float pitch;
	float yaw;
}angel_t;

typedef struct 
{
  float q0;
  float q1;
  float q2;
  float q3;
} Quaternion;
/* Return the difference of two vectors, (vector1 - vector2). */
vec3d vdiff(const vec3d vector1, const vec3d vector2);

/* Return the sum of two vectors. */
vec3d vsum(const vec3d vector1, const vec3d vector2);

/* Multiply vector by a number. */
vec3d vmul(const vec3d vector, const double n);

/* Divide vector by a number. */
vec3d vdiv(const vec3d vector, const double n);

/* Return the Euclidean norm. */
double vdist(const vec3d v1, const vec3d v2);

/* Return the Euclidean norm. */
double vnorm(const vec3d vector);

/* Return the dot product of two vectors. */
double dot(const vec3d vector1, const vec3d vector2);

/* Replace vector with its cross product with another vector. */
vec3d cross(const vec3d vector1, const vec3d vector2);

int GetLocation(vec3d *best_solution, vec3d* anchorArray, int *distanceArray);

double vdist(const vec3d v1, const vec3d v2);
void IMUupdate(double gx, double gy, double gz, double ax, double ay, double az);
//加速度，角速度获取四元数
Quaternion GetAngle(const ImuData_t *pMpu, float dt);
#endif
