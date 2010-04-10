#include <math.h>
#include "main.h"

/* Needs "-std=c99" to compile. */

static real_t dot(const real_t *v, const real_t *u, uint_t dim)
{
    real_t res = 0;
    for (uint_t k = 0; k < dim; ++k)
        res = res + v[k] * u[k];
    return res;
}

static real_t normalize(real_t *v, uint_t dim)
{
    real_t norm = sqrt(dot(v, v, dim));

    real_t inv_norm = 1.0 / norm;
    for (uint_t k = 0; k < dim; ++k)
        v[k] = v[k] * inv_norm;

    return norm;
}

static void mpe(real_t *R, real_t *gamma, uint_t num)
{
    /* TODO */
}

static void rre(real_t *R, real_t *gamma, uint_t num)
{
    /* TODO */
}

static void mgs(real_t *vectors, real_t *R, uint_t dim, uint_t num)
{
    for (uint_t k = 0; k < num * num; ++k)
        R[k] = 0;

    for (uint_t i = 0; i < num; ++i) {
        real_t *v = vectors + i * dim;

        for (uint_t j = 0; j < i; ++j) {

            real_t *u = vectors + j * dim;
            real_t r = dot(v, u, dim);

            for (uint_t k = 0; k < dim; ++k)
                v[k] = v[k] - r * u[k];

            R[j*dim + i] = r;
        }
        R[i*dim + i] = normalize(v, dim);
    }
}

int DLL_EXPORT extrapolate(real_t *vectors, uint_t dim, uint_t num, enum EXTR_METHOD method)
{
    real_t *const v0 = vectors;

    /* Compute differences */
    for (uint_t i = num - 1; i > 0; --i) {
        real_t *v = vectors + i * dim;
        real_t *u = v - dim;
        for (uint_t k = 0; k < dim; ++k)
            v[k] = v[k] - u[k];
    }
    vectors = vectors + dim;
    num = num - 1;

    /* QR decomposition */
    real_t R[num * num];
    mgs(vectors, R, dim, num);

    /* Compute coefficients */
    real_t gamma[num];
    switch (method) {
        case MPE:	mpe(R, gamma, num); break;
        case RRE:  	rre(R, gamma, num); break;
        default:    return -1;
    }

    /* xi = 1 - cumsum(gamma) */
    for (uint_t i = 1; i < num; ++i)
        gamma[i] = gamma[i] + gamma[i-1];
    real_t xi[num];
    for (uint_t i = 0; i < num; ++i)
        xi[i] = 1 - gamma[i];

    /* eta = R * xi */
    real_t eta[num];
    for (uint_t i = 0; i < num; ++i)
        eta[i] = dot(R + i*num, xi, num);

    /* Compute the solution: v0 + Q * eta */
    for (uint_t i = 0; i < num; ++i) {
        real_t *q = vectors + i * dim;
        for (uint_t k = 0; k < dim; ++k) {
            v0[k] = v0[k] + q[k] * eta[i];
        }
    }
    return 0;
}
