// -------------------------------------------------------------------------------------------------------------------
//
//  File: trilateration.cpp
//
//  Copyright 2016 (c) Decawave Ltd, Dublin, Ireland.
//
//  All rights reserved.
//
//  Author:
//
// -------------------------------------------------------------------------------------------------------------------

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "Uwb_Location/trilateration.h"
#include "ceres/ceres.h"
#include "iostream"

using namespace std;
vec3d prev_res;
double cost1 = 1;

/* Using Ceres to perform optimation */
// struct TRILATERATION_COST {
//     TRILATERATION_COST(vec3d p1, vec3d p2, vec3d p3, vec3d p4, double r1, double r2, double r3, double r4) : p1_(p1), p2_(p2), p3_(p3), p4_(p4), r1_(r1), r2_(r2), r3_(r3), r4_(r4) {};

//     template <typename T>
//     bool operator()(const T* const xyz, T* residual) const
//     {
//         residual[0] = (xyz[0] - T(p1_.x)) * (xyz[0] - T(p1_.x)) + (xyz[1] - T(p1_.y)) * (xyz[1] - T(p1_.y)) + (xyz[2] - T(p1_.z)) * (xyz[2] - T(p1_.z)) - T(r1_) * T(r1_);
//         residual[1] = (xyz[0] - T(p2_.x)) * (xyz[0] - T(p2_.x)) + (xyz[1] - T(p2_.y)) * (xyz[1] - T(p2_.y)) + (xyz[2] - T(p2_.z)) * (xyz[2] - T(p2_.z)) - T(r2_) * T(r2_);
//         residual[2] = (xyz[0] - T(p3_.x)) * (xyz[0] - T(p3_.x)) + (xyz[1] - T(p3_.y)) * (xyz[1] - T(p3_.y)) + (xyz[2] - T(p3_.z)) * (xyz[2] - T(p3_.z)) - T(r3_) * T(r3_);
//         residual[3] = (xyz[0] - T(p4_.x)) * (xyz[0] - T(p4_.x)) + (xyz[1] - T(p4_.y)) * (xyz[1] - T(p4_.y)) + (xyz[2] - T(p4_.z)) * (xyz[2] - T(p4_.z)) - T(r4_) * T(r4_);
//         return true;
//     }

//     const double r1_, r2_, r3_, r4_;
//     const vec3d p1_, p2_, p3_, p4_;
// };

struct TRILATERATION_COST {
    TRILATERATION_COST(vec3d p, double r)
        : p_(p)
        , r_(r) {};

    template <typename T>
    bool operator()(const T* const xyz, T* residual) const
    {
        residual[0] = sqrt((xyz[0] - T(p_.x)) * (xyz[0] - T(p_.x)) + (xyz[1] - T(p_.y)) * (xyz[1] - T(p_.y)) + (xyz[2] - T(p_.z)) * (xyz[2] - T(p_.z))) - T(r_);
        return true;
    }

    const double r_;
    const vec3d p_;
};

/* Largest nonnegative number still considered zero */
#define MAXZERO 0.001

#define ERR_TRIL_CONCENTRIC -1
#define ERR_TRIL_COLINEAR_2SOLUTIONS -2
#define ERR_TRIL_SQRTNEGNUMB -3
#define ERR_TRIL_NOINTERSECTION_SPHERE4 -4
#define ERR_TRIL_NEEDMORESPHERE -5
#define ERR_TRIL_CANNOTDETERMINE -6

#define CM_ERR_ADDED (10) // was 5

/* Return the difference of two vectors, (vector1 - vector2). */
vec3d vdiff(const vec3d vector1, const vec3d vector2)
{
    vec3d v;
    v.x = vector1.x - vector2.x;
    v.y = vector1.y - vector2.y;
    v.z = vector1.z - vector2.z;
    return v;
}

/* Return the sum of two vectors. */
vec3d vsum(const vec3d vector1, const vec3d vector2)
{
    vec3d v;
    v.x = vector1.x + vector2.x;
    v.y = vector1.y + vector2.y;
    v.z = vector1.z + vector2.z;
    return v;
}

/* Multiply vector by a number. */
vec3d vmul(const vec3d vector, const double n)
{
    vec3d v;
    v.x = vector.x * n;
    v.y = vector.y * n;
    v.z = vector.z * n;
    return v;
}

/* Divide vector by a number. */
vec3d vdiv(const vec3d vector, const double n)
{
    vec3d v;
    v.x = vector.x / n;
    v.y = vector.y / n;
    v.z = vector.z / n;
    return v;
}

/* Return the Euclidean norm. */
double vdist(const vec3d v1, const vec3d v2)
{
    double xd = v1.x - v2.x;
    double yd = v1.y - v2.y;
    double zd = v1.z - v2.z;
    return sqrt(xd * xd + yd * yd + zd * zd);
}

/* Return the Euclidean norm. */
double vnorm(const vec3d vector)
{
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

/* Return the dot product of two vectors. */
double dot(const vec3d vector1, const vec3d vector2)
{
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

/* Replace vector with its cross product with another vector. */
vec3d cross(const vec3d vector1, const vec3d vector2)
{
    vec3d v;
    v.x = vector1.y * vector2.z - vector1.z * vector2.y;
    v.y = vector1.z * vector2.x - vector1.x * vector2.z;
    v.z = vector1.x * vector2.y - vector1.y * vector2.x;
    return v;
}

/* Return the GDOP (Geometric Dilution of Precision) rate between 0-1.
 * Lower GDOP rate means better precision of intersection.
 */
double gdoprate(const vec3d tag, const vec3d p1, const vec3d p2, const vec3d p3)
{
    vec3d ex, t1, t2, t3;
    double h, gdop1, gdop2, gdop3, result;

    ex = vdiff(p1, tag);
    h = vnorm(ex);
    t1 = vdiv(ex, h);

    ex = vdiff(p2, tag);
    h = vnorm(ex);
    t2 = vdiv(ex, h);

    ex = vdiff(p3, tag);
    h = vnorm(ex);
    t3 = vdiv(ex, h);

    gdop1 = fabs(dot(t1, t2));
    gdop2 = fabs(dot(t2, t3));
    gdop3 = fabs(dot(t3, t1));

    if (gdop1 < gdop2)
        result = gdop2;
    else
        result = gdop1;
    if (result < gdop3)
        result = gdop3;

    return result;
}

/* Intersecting a sphere sc with radius of r, with a line p1-p2.
 * Return zero if successful, negative error otherwise.
 * mu1 & mu2 are constant to find points of intersection.
 */
int sphereline(const vec3d p1, const vec3d p2, const vec3d sc, double r, double* const mu1, double* const mu2)
{
    double a, b, c;
    double bb4ac;
    vec3d dp;

    dp.x = p2.x - p1.x;
    dp.y = p2.y - p1.y;
    dp.z = p2.z - p1.z;

    a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

    b = 2 * (dp.x * (p1.x - sc.x) + dp.y * (p1.y - sc.y) + dp.z * (p1.z - sc.z));

    c = sc.x * sc.x + sc.y * sc.y + sc.z * sc.z;
    c += p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;
    c -= 2 * (sc.x * p1.x + sc.y * p1.y + sc.z * p1.z);
    c -= r * r;

    bb4ac = b * b - 4 * a * c;

    if (fabs(a) == 0 || bb4ac < 0) {
        *mu1 = 0;
        *mu2 = 0;
        return -1;
    }

    *mu1 = (-b + sqrt(bb4ac)) / (2 * a);
    *mu2 = (-b - sqrt(bb4ac)) / (2 * a);

    return 0;
}

/* Return TRIL_3SPHERES if it is performed using 3 spheres and return
 * TRIL_4SPHERES if it is performed using 4 spheres
 * For TRIL_3SPHERES, there are two solutions: result1 and result2
 * For TRIL_4SPHERES, there is only one solution: best_solution
 *
 * Return negative number for other errors
 *
 * To force the function to work with only 3 spheres, provide a duplicate of
 * any sphere at any place among p1, p2, p3 or p4.
 *
 * The last parameter is the largest nonnegative number considered zero;
 * it is somewhat analogous to machine epsilon (but inclusive).
 */
int trilateration(vec3d* const result1,
    vec3d* const result2,
    vec3d* const best_solution,
    const vec3d p1, const double r1,
    const vec3d p2, const double r2,
    const vec3d p3, const double r3,
    const vec3d p4, const double r4,
    const double maxzero)
{
    vec3d ex, ey, ez, t1, t2, t3;
    double h, i, j, x, y, z, t;
    double mu1, mu2, mu;
    int result;

    /*********** FINDING TWO POINTS FROM THE FIRST THREE SPHERES **********/

    // if there are at least 2 concentric spheres within the first 3 spheres
    // then the calculation may not continue, drop it with error -1

    /* h = |p3 - p1|, ex = (p3 - p1) / |p3 - p1| */
    ex = vdiff(p3, p1); // vector p13
    h = vnorm(ex); // scalar p13
    if (h <= maxzero) {
        /* p1 and p3 are concentric, not good to obtain a precise intersection point */
        // printf("concentric13 return -1\n");
        return ERR_TRIL_CONCENTRIC;
    }

    /* h = |p3 - p2|, ex = (p3 - p2) / |p3 - p2| */
    ex = vdiff(p3, p2); // vector p23
    h = vnorm(ex); // scalar p23
    if (h <= maxzero) {
        /* p2 and p3 are concentric, not good to obtain a precise intersection point */
        // printf("concentric23 return -1\n");
        return ERR_TRIL_CONCENTRIC;
    }

    /* h = |p2 - p1|, ex = (p2 - p1) / |p2 - p1| */
    ex = vdiff(p2, p1); // vector p12
    h = vnorm(ex); // scalar p12
    if (h <= maxzero) {
        /* p1 and p2 are concentric, not good to obtain a precise intersection point */
        // printf("concentric12 return -1\n");
        return ERR_TRIL_CONCENTRIC;
    }
    ex = vdiv(ex, h); // unit vector ex with respect to p1 (new coordinate system)

    /* t1 = p3 - p1, t2 = ex (ex . (p3 - p1)) */
    t1 = vdiff(p3, p1); // vector p13
    i = dot(ex, t1); // the scalar of t1 on the ex direction
    t2 = vmul(ex, i); // colinear vector to p13 with the length of i

    /* ey = (t1 - t2), t = |t1 - t2| */
    ey = vdiff(t1, t2); // vector t21 perpendicular to t1
    t = vnorm(ey); // scalar t21
    if (t > maxzero) {
        /* ey = (t1 - t2) / |t1 - t2| */
        ey = vdiv(ey, t); // unit vector ey with respect to p1 (new coordinate system)

        /* j = ey . (p3 - p1) */
        j = dot(ey, t1); // scalar t1 on the ey direction
    } else
        j = 0.0;

    /* Note: t <= maxzero implies j = 0.0. */
    if (fabs(j) <= maxzero) {

        /* Is point p1 + (r1 along the axis) the intersection? */
        t2 = vsum(p1, vmul(ex, r1));
        if (fabs(vnorm(vdiff(p2, t2)) - r2) <= maxzero && fabs(vnorm(vdiff(p3, t2)) - r3) <= maxzero) {
            /* Yes, t2 is the only intersection point. */
            if (result1)
                *result1 = t2;
            if (result2)
                *result2 = t2;
            return TRIL_3SPHERES;
        }

        /* Is point p1 - (r1 along the axis) the intersection? */
        t2 = vsum(p1, vmul(ex, -r1));
        if (fabs(vnorm(vdiff(p2, t2)) - r2) <= maxzero && fabs(vnorm(vdiff(p3, t2)) - r3) <= maxzero) {
            /* Yes, t2 is the only intersection point. */
            if (result1)
                *result1 = t2;
            if (result2)
                *result2 = t2;
            return TRIL_3SPHERES;
        }
        /* p1, p2 and p3 are colinear with more than one solution */
        return ERR_TRIL_COLINEAR_2SOLUTIONS;
    }

    /* ez = ex x ey */
    ez = cross(ex, ey); // unit vector ez with respect to p1 (new coordinate system)

    x = (r1 * r1 - r2 * r2) / (2 * h) + h / 2;
    y = (r1 * r1 - r3 * r3 + i * i) / (2 * j) + j / 2 - x * i / j;
    z = r1 * r1 - x * x - y * y;
    std::cout << "r1 = " << r1 << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    if (z < -maxzero) {
        /* The solution is invalid, square root of negative number */
        return ERR_TRIL_SQRTNEGNUMB;
    } else if (z > 0.0)
        z = sqrt(z);
    else
        z = 0.0;

    /* t2 = p1 + x ex + y ey */
    t2 = vsum(p1, vmul(ex, x));
    t2 = vsum(t2, vmul(ey, y));

    /* result1 = p1 + x ex + y ey + z ez */
    if (result1)
        *result1 = vsum(t2, vmul(ez, z));

    /* result1 = p1 + x ex + y ey - z ez */
    if (result2)
        *result2 = vsum(t2, vmul(ez, -z));

    /*********** END OF FINDING TWO POINTS FROM THE FIRST THREE SPHERES **********/
    /********* RESULT1 AND RESULT2 ARE SOLUTIONS, OTHERWISE RETURN ERROR *********/

    /************* FINDING ONE SOLUTION BY INTRODUCING ONE MORE SPHERE ***********/

    // check for concentricness of sphere 4 to sphere 1, 2 and 3
    // if it is concentric to one of them, then sphere 4 cannot be used
    // to determine the best solution and return -1

    /* h = |p4 - p1|, ex = (p4 - p1) / |p4 - p1| */
    ex = vdiff(p4, p1); // vector p14
    h = vnorm(ex); // scalar p14
    if (h <= maxzero) {
        /* p1 and p4 are concentric, not good to obtain a precise intersection point */
        //     // printf("concentric14 return 0\n");
        return TRIL_3SPHERES;
    }
    /* h = |p4 - p2|, ex = (p4 - p2) / |p4 - p2| */
    ex = vdiff(p4, p2); // vector p24
    h = vnorm(ex); // scalar p24
    if (h <= maxzero) {
        /* p2 and p4 are concentric, not good to obtain a precise intersection point */
        // printf("concentric24 return 0\n");
        return TRIL_3SPHERES;
    }
    /* h = |p4 - p3|, ex = (p4 - p3) / |p4 - p3| */
    ex = vdiff(p4, p3); // vector p34
    h = vnorm(ex); // scalar p34
    if (h <= maxzero) {
        /* p3 and p4 are concentric, not good to obtain a precise intersection point */
        // printf("concentric34 return 0\n");
        return TRIL_3SPHERES;
    }

    /* With Ceres */
    // double xyz[3] = {result1->x, result1->y, result1->z};
    // ceres::Problem problem;
    // problem.AddResidualBlock(new ceres::AutoDiffCostFunction<TRILATERATION_COST, 4, 3>(new TRILATERATION_COST(p1, p2, p3, p4, r1, r2, r3, r4)), new ceres::HuberLoss(0.5), xyz);
    // problem.SetParameterLowerBound(xyz, 0, 0.3);
    // problem.SetParameterLowerBound(xyz, 1, 0.3);
    // problem.SetParameterLowerBound(xyz, 2, 0.3);
    // ceres::Solver::Options options;
    // options.linear_solver_type=ceres::DENSE_QR;
    // options.minimizer_progress_to_stdout=true;
    // ceres::Solver::Summary summary;
    // ceres::Solve(options, &problem, &summary);
    // best_solution->x = xyz[0];
    // best_solution->y = xyz[1];
    // best_solution->z = xyz[2];
    // return TRIL_4SPHERES;

    // std::cout << "result2 " << result2->x << " " << result2->y << " " << result2->z << std::endl;
    double xyz[3];
    vec3d tmp;

    // if ((gdoprate(*result1, p1, p2, p4) < gdoprate(*result2, p1, p2, p4)) &&
    //  (gdoprate(*result2, p1, p2, p4) < gdoprate(prev_res, p1, p2, p4)+0.1))  {
    //     tmp = vdiv(vsum(*result1, prev_res), 2);
    // } 
    // else if ((gdoprate(*result2, p1, p2, p4) < gdoprate(*result1, p1, p2, p4)) &&
    //  (gdoprate(*result1, p1, p2, p4) < gdoprate(prev_res, p1, p2, p4)+0.1))  {
    //     tmp = vdiv(vsum(*result2, prev_res), 2);
    // } 
    // // one of them is larger than prev_res
    // else if (gdoprate(*result2, p1, p2, p4) < gdoprate(prev_res, p1, p2, p4)+0.025) {
    //     tmp = vdiv(vsum(*result2, prev_res), 2);
    // } 
    // else if (gdoprate(*result1, p1, p2, p4) < gdoprate(prev_res, p1, p2, p4)+0.025) {
    //     tmp = vdiv(vsum(*result1, prev_res), 2);
    // } 
    // // both of them are larger than prev_res
    // else {
    //     tmp = prev_res;
    // }

    double compare1, compare2;
    compare1 = vnorm(vdiff(*result1, p4));
    compare2 = vnorm(vdiff(*result2, p4));
    if(fabs(compare1-compare2) <= 0.2){
        if(prev_res.x == 10000.0){
            tmp = vdiv(vsum(*result1, *result2), 2);
        }
        else{
            tmp = prev_res;
        }
    }

    std::cout<< fabs(compare1 - r4) << " " << fabs(compare2 - r4)  << std::endl;
    if(fabs(compare1 - r4) < fabs(compare2 - r4)){
        tmp = vdiv(vsum(*result1, prev_res), 2);
        std::cout<< "error: " << fabs(compare1 - r4) << " " << fabs(compare2 - r4)  << std::endl;
        std::cout<< "choose result1"  << std::endl;
    }
    else{
        tmp = vdiv(vsum(*result2, prev_res), 2);
        std::cout<< "choose result2"  << std::endl;
    }
    std::cout<< "result1 = " << result1->x << " " << result1->y << " " << result1->z << std::endl;
    std::cout<< "result2 = " << result2->x << " " << result2->y << " " << result2->z << std::endl;
    std::cout<< "the primal point at(" << tmp.x << ", " << tmp.y << ", " << tmp.z << " )" << std::endl;
    xyz[0] = tmp.x;
    xyz[1] = tmp.y;
    xyz[2] = tmp.z;
    ceres::Problem problem;
    vec3d p[4] = { p1, p2, p3, p4 };
    double r[4] = { r1, r2, r3, r4 };
    for (int i = 0; i < 4; ++i) {
        problem.AddResidualBlock(new ceres::AutoDiffCostFunction<TRILATERATION_COST, 1, 3>(new TRILATERATION_COST(p[i], r[i])), nullptr, xyz);
        //problem.SetParameterLowerBound(xyz, 2, 0.);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    options.max_num_iterations = 200;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cost1 = summary.final_cost;
    std::cout << "Total cost of this iteration: " << summary.final_cost << std::endl;
    best_solution->x = xyz[0];
    best_solution->y = xyz[1];
    best_solution->z = xyz[2];
    prev_res = *best_solution;
    std::cout <<"the optimal point at(" << best_solution->x << ", " << best_solution->y << ", " << best_solution->z << " )" << std::endl;
    return TRIL_4SPHERES;

    /******** END OF FINDING ONE SOLUTION BY INTRODUCING ONE MORE SPHERE *********/
}

/* This function calls trilateration to get the best solution.
 *
 * If any three spheres does not produce valid solution,
 * then each distance is increased to ensure intersection to happens.
 *
 * Return the selected trilateration mode between TRIL_3SPHERES or TRIL_4SPHERES
 * For TRIL_3SPHERES, there are two solutions: solution1 and solution2
 * For TRIL_4SPHERES, there is only one solution: best_solution
 *
 * nosolution_count = the number of failed attempt before intersection is found
 * by increasing the sphere diameter.
 */
int deca_3dlocate(vec3d* const solution1,
    vec3d* const solution2,
    vec3d* const best_solution,
    int* const nosolution_count,
    double* const best_3derror,
    double* const best_gdoprate,
    vec3d p1, double r1,
    vec3d p2, double r2,
    vec3d p3, double r3,
    vec3d p4, double r4,
    int* combination)
{
    vec3d o1, o2, solution, ptemp;
    // vec3d    solution_compare1, solution_compare2;
    double /*error_3dcompare1, error_3dcompare2,*/ rtemp;
    double best_cost = 1, cost2 = 1;
    double gdoprate_compare1, gdoprate_compare2;
    double ovr_r1, ovr_r2, ovr_r3, ovr_r4;
    int overlook_count, combination_counter;
    int trilateration_errcounter, trilateration_mode34;
    int success, concentric, result;

    trilateration_errcounter = 0;
    trilateration_mode34 = 0;

    combination_counter = 4; /* four spheres combination */

    *best_gdoprate = 1; /* put the worst gdoprate init */
    gdoprate_compare1 = 1;
    gdoprate_compare2 = 1;
    // solution_compare1.x = 0; solution_compare1.y = 0; solution_compare1.z = 0;
    // error_3dcompare1 = 0;

    do {
        success = 0;
        concentric = 0;
        overlook_count = 0;
        ovr_r1 = r1;
        ovr_r2 = r2;
        ovr_r3 = r3;
        ovr_r4 = r4;

        do {

            result = trilateration(&o1, &o2, &solution, p1, ovr_r1, p2, ovr_r2, p3, ovr_r3, p4, ovr_r4, MAXZERO);

            switch (result) {
            case TRIL_3SPHERES: // 3 spheres are used to get the result
                trilateration_mode34 = TRIL_3SPHERES;
                success = 1;
                break;

            case TRIL_4SPHERES: // 4 spheres are used to get the result
                trilateration_mode34 = TRIL_4SPHERES;
                success = 1;
                break;

            case ERR_TRIL_CONCENTRIC:
                concentric = 1;
                break;

            default: // any other return value goes here
                ovr_r1 += 0.10;
                ovr_r2 += 0.10;
                ovr_r3 += 0.10;
                ovr_r4 += 0.10;
                overlook_count++;
                break;
            }

        } while (!success && (overlook_count <= CM_ERR_ADDED) && !concentric);

        if (success) {
            switch (result) {
            case TRIL_3SPHERES:
                *solution1 = o1;
                *solution2 = o2;
                *nosolution_count = overlook_count;

                combination_counter = 0;
                break;

            case TRIL_4SPHERES:
                /* calculate the new gdop */
                std::cout <<"cost1 = " << cost1 << ", " << "cost2 = " << cost2 << std::endl;

                if (cost1 <= cost2) {
                    *solution1 = o1;
                    *solution2 = o2;
                    *best_solution = solution;
                    *best_3derror = sqrt((vnorm(vdiff(solution, p1)) - r1) * (vnorm(vdiff(solution, p1)) - r1) 
                                    + (vnorm(vdiff(solution, p2)) - r2) * (vnorm(vdiff(solution, p2)) - r2) 
                                    + (vnorm(vdiff(solution, p3)) - r3) * (vnorm(vdiff(solution, p3)) - r3) 
                                    + (vnorm(vdiff(solution, p4)) - r4) * (vnorm(vdiff(solution, p4)) - r4));
                    best_cost = cost1;
                    cost2 = cost1;
                }
                *combination = 5 - combination_counter;

                ptemp = p1;
                p1 = p2;
                p2 = p3;
                p3 = p4;
                p4 = ptemp;
                rtemp = r1;
                r1 = r2;
                r2 = r3;
                r3 = r4;
                r4 = rtemp;
                combination_counter--;
                break;

            default:
                break;
            }
        } else // if (success)
        {
            trilateration_errcounter++;
            combination_counter--;
        }

    } while (combination_counter);

    // printf("trilateration_errcounter = %d\n", trilateration_errcounter);
    //  if it gives error for all 4 sphere combinations then no valid result is given
    //  otherwise return the trilateration mode used
    std::cout <<"the optimal point at(" << solution.x << ", " << solution.y << ", " << solution.z << " )" << std::endl;

    if (trilateration_errcounter >= 4)
        return result;
    else
        return trilateration_mode34;
}

struct num {
    int anc_ID;
    int distance;
} valid_anc_num[(MAX_AHCHOR_NUMBER + 1)];

int cmp(const void* m, const void* n) // 定义返回值返回方式
{
    return ((struct num*)m)->distance - ((struct num*)n)->distance;
}

int GetLocation(vec3d* best_solution, vec3d* anchorArray, int* distanceArray)
{

    vec3d o1, o2, p1, p2, p3, p4;
    double r1 = 0, r2 = 0, r3 = 0, r4 = 0, best_3derror, best_gdoprate;
    int result;
    int error, combination;
    int valid_anc_count = 0;
    int j = 0;
    int use3anc = 0;

    for (int i = 0; i < (MAX_AHCHOR_NUMBER + 1); i++) // 清空结构体数组
    {
        valid_anc_num[i].anc_ID = 0;
        valid_anc_num[i].distance = 0;
    }

    for (int i = 0; i < MAX_AHCHOR_NUMBER; i++) // 验证几个有效距离值
    {
        if (distanceArray[i] > 0) {
            valid_anc_count++;
            valid_anc_num[j].anc_ID = i; // 记录有效基站编号
            valid_anc_num[j].distance = distanceArray[i]; // 记录有效基站距离
            j++;
        }
    }

    if (valid_anc_count < 3) {
        return -1;
        puts("err1");
    }

    else if (valid_anc_count == 3) // 直接执行三基站定位
    {
        use3anc = 1;
        /* Anchors coordinate */
        p1.x = anchorArray[valid_anc_num[0].anc_ID].x;
        p1.y = anchorArray[valid_anc_num[0].anc_ID].y;
        p1.z = anchorArray[valid_anc_num[0].anc_ID].z;
        p2.x = anchorArray[valid_anc_num[1].anc_ID].x;
        p2.y = anchorArray[valid_anc_num[1].anc_ID].y;
        p2.z = anchorArray[valid_anc_num[1].anc_ID].z;
        p3.x = anchorArray[valid_anc_num[2].anc_ID].x;
        p3.y = anchorArray[valid_anc_num[2].anc_ID].y;
        p3.z = anchorArray[valid_anc_num[2].anc_ID].z;
        p4.x = p1.x;
        p4.y = p1.y;
        p4.z = p1.z;

        r1 = (double)distanceArray[valid_anc_num[0].anc_ID] / 1000.0;
        r2 = (double)distanceArray[valid_anc_num[1].anc_ID] / 1000.0;
        r3 = (double)distanceArray[valid_anc_num[2].anc_ID] / 1000.0;
        r4 = r1;

        // printf("anc1=%d,anc2=%d,anc3=%d\n",valid_anc_num[0],valid_anc_num[1],valid_anc_num[2]);
        // printf("dis1=%f,dis2=%f,dis3=%f\n",r1,r2,r3);
        // printf("P1X=%f,P1Y=%f,P1Z=%f\n",p1.x,p1.y,p1.z);
    }

    else if (valid_anc_count == 4) // 直接执行4基站定位
    {
        /* Anchors coordinate */
        p1.x = anchorArray[valid_anc_num[0].anc_ID].x;
        p1.y = anchorArray[valid_anc_num[0].anc_ID].y;
        p1.z = anchorArray[valid_anc_num[0].anc_ID].z;
        p2.x = anchorArray[valid_anc_num[1].anc_ID].x;
        p2.y = anchorArray[valid_anc_num[1].anc_ID].y;
        p2.z = anchorArray[valid_anc_num[1].anc_ID].z;
        p3.x = anchorArray[valid_anc_num[2].anc_ID].x;
        p3.y = anchorArray[valid_anc_num[2].anc_ID].y;
        p3.z = anchorArray[valid_anc_num[2].anc_ID].z;
        p4.x = anchorArray[valid_anc_num[3].anc_ID].x;
        p4.y = anchorArray[valid_anc_num[3].anc_ID].y;
        p4.z = anchorArray[valid_anc_num[3].anc_ID].z;

        r1 = (double)distanceArray[valid_anc_num[0].anc_ID] / 1000.0;
        r2 = (double)distanceArray[valid_anc_num[1].anc_ID] / 1000.0;
        r3 = (double)distanceArray[valid_anc_num[2].anc_ID] / 1000.0;
        r4 = (double)distanceArray[valid_anc_num[3].anc_ID] / 1000.0;

    }

    // valid_anc_count 有效基站个数
    // valid_anc_num[0] 有效基站编号
    else if (valid_anc_count > 4) // 执行基站选取机制，选取最近的4个基站1234进行计算
    {
        qsort(valid_anc_num, (valid_anc_count + 1), sizeof(valid_anc_num[0]), cmp); // 将有效距离值进行从小到大排序
        for (int i = 1; i <= valid_anc_count; i++) // 输出结果
            printf("No%d DIS=%d,ID=A%d\n", i, valid_anc_num[i].distance, valid_anc_num[i].anc_ID);

        p1.x = anchorArray[valid_anc_num[1].anc_ID].x;
        p1.y = anchorArray[valid_anc_num[1].anc_ID].y;
        p1.z = anchorArray[valid_anc_num[1].anc_ID].z;
        p2.x = anchorArray[valid_anc_num[2].anc_ID].x;
        p2.y = anchorArray[valid_anc_num[2].anc_ID].y;
        p2.z = anchorArray[valid_anc_num[2].anc_ID].z;
        p3.x = anchorArray[valid_anc_num[3].anc_ID].x;
        p3.y = anchorArray[valid_anc_num[3].anc_ID].y;
        p3.z = anchorArray[valid_anc_num[3].anc_ID].z;
        p4.x = anchorArray[valid_anc_num[4].anc_ID].x;
        p4.y = anchorArray[valid_anc_num[4].anc_ID].y;
        p4.z = anchorArray[valid_anc_num[4].anc_ID].z;

        r1 = (double)distanceArray[valid_anc_num[1].anc_ID] / 1000.0;
        r2 = (double)distanceArray[valid_anc_num[2].anc_ID] / 1000.0;
        r3 = (double)distanceArray[valid_anc_num[3].anc_ID] / 1000.0;
        r4 = (double)distanceArray[valid_anc_num[4].anc_ID] / 1000.0;

        printf("use1=A%d,use2=A%d,use3=A%d,use4=A%d\n", valid_anc_num[1].anc_ID, valid_anc_num[2].anc_ID, valid_anc_num[3].anc_ID, valid_anc_num[4].anc_ID);
    }

    result = deca_3dlocate(&o1, &o2, best_solution, &error, &best_3derror, &best_gdoprate,
        p1, r1, p2, r2, p3, r3, p4, r4, &combination);

    if ((result == 0) && (valid_anc_count > 4)) // 多于4基站选取后计算失败，把第1舍掉用第2345计算
    {
        puts("Second calculation");
        p1.x = anchorArray[valid_anc_num[2].anc_ID].x;
        p1.y = anchorArray[valid_anc_num[2].anc_ID].y;
        p1.z = anchorArray[valid_anc_num[2].anc_ID].z;
        p2.x = anchorArray[valid_anc_num[3].anc_ID].x;
        p2.y = anchorArray[valid_anc_num[3].anc_ID].y;
        p2.z = anchorArray[valid_anc_num[3].anc_ID].z;
        p3.x = anchorArray[valid_anc_num[4].anc_ID].x;
        p3.y = anchorArray[valid_anc_num[4].anc_ID].y;
        p3.z = anchorArray[valid_anc_num[4].anc_ID].z;
        p4.x = anchorArray[valid_anc_num[5].anc_ID].x;
        p4.y = anchorArray[valid_anc_num[5].anc_ID].y;
        p4.z = anchorArray[valid_anc_num[5].anc_ID].z;

        r1 = (double)distanceArray[valid_anc_num[2].anc_ID] / 1000.0;
        r2 = (double)distanceArray[valid_anc_num[3].anc_ID] / 1000.0;
        r3 = (double)distanceArray[valid_anc_num[4].anc_ID] / 1000.0;
        r4 = (double)distanceArray[valid_anc_num[5].anc_ID] / 1000.0;

        printf("use1=A%d,use2=A%d,use3=A%d,use4=A%d\n", valid_anc_num[2].anc_ID, valid_anc_num[3].anc_ID, valid_anc_num[4].anc_ID, valid_anc_num[5].anc_ID);

        result = deca_3dlocate(&o1, &o2, best_solution, &error, &best_3derror, &best_gdoprate,
            p1, r1, p2, r2, p3, r3, p4, r4, &combination);
    }
    

    if (result >= 0) {
        return result;
    }

    // if (result >= 0) {
    //     if (o1.z <= o2.z)
    //         best_solution->z = o1.z;
    //     else
    //         best_solution->z = o2.z;
    //     if (use3anc == 1 || result == TRIL_3SPHERES) {
    //         if (o1.z < p1.z)
    //             *best_solution = o1;
    //         else
    //             *best_solution = o2; // assume tag is below the anchors (1, 2, and 3)
    //     }

    //     return result;
    // }
    return result;
}
