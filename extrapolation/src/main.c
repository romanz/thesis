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

static real_t sum(real_t *v, uint_t dim)
{
    real_t result = 0;

    for (uint_t k = 0; k < dim; ++k)
        result = result + v[k];

    return result;
}

#define IDX(I, J, DIM) ((I) + (J)*(DIM))

void DLL_EXPORT solver(real_t *sol, const real_t *A, real_t *b, uint_t dim, bool_t upper_diagonal)
{
    if (upper_diagonal) {
        for (uint_t i = dim; i > 0;) {
            --i;
            real_t x = b[i] / A[IDX(i, i, dim)];
            for (uint_t k = 0; k < i; ++k)
                b[k] = b[k] - x * A[IDX(k, i, dim)];
            sol[i] = x;
        }
    } else {
        for (uint_t i = 0; i < dim;) {
            real_t x = b[i] / A[IDX(i, i, dim)];
            for (uint_t k = i; k < dim; ++k)
                b[k] = b[k] - x * A[IDX(k, i, dim)];
            sol[i] = x;
            ++i;
        }
    }
}

static real_t mpe(real_t *R, real_t *gamma, uint_t order)
{
    real_t A[(order - 1)*(order - 1)];
    for (uint_t i = 0; i < order - 1; ++i)
        for (uint_t j = 0; j < order - 1; ++j)
            A[IDX(i, j, order - 1)] = R[IDX(i, j, order)];

    real_t b[order];
    for (uint_t i = 0; i < order - 1; ++i)
        b[i] = -R[IDX(i, order - 1, order)];

    solver(b, A, b, order - 1, 1);
    b[order - 1] = 1;

    real_t s = 1.0 / sum(b, order);
    for (uint_t k = 0; k < order; ++k)
        gamma[k] = b[k] * s;
    return fabs(gamma[order - 1]) * R[IDX(order - 1, order - 1, order)];
}

static void transpose(real_t *A, uint_t dim)
{
    for (uint_t i = 0; i < dim; ++i) {
        for (uint_t j = 0; j < i; ++j) {
            real_t x = A[IDX(i, j, dim)];
            A[IDX(i, j, dim)] = A[IDX(j, i, dim)];
            A[IDX(j, i, dim)] = x;
        }
    }
}

static real_t rre(real_t *R, real_t *gamma, uint_t order)
{
    real_t b[order];
    for (uint_t k = 0; k < order; ++k)
        b[k] = 1;

    transpose(R, order);
    solver(b, R, b, order, 1);

    transpose(R, order);
    solver(b, R, b, order, 0);

    real_t lambda = 0;
    for (uint_t k = 0; k < order; ++k)
        lambda = lambda + b[k];
    lambda = 1.0 / lambda;
    for (uint_t k = 0; k < order; ++k)
        gamma[k] = b[k] * lambda;
    return sqrt(lambda);
}

static void mgs(real_t *vectors, real_t *R, uint_t dim, uint_t order)
{
    for (uint_t k = 0; k < order * order; ++k)
        R[k] = 0;

    for (uint_t i = 0; i < order; ++i) {
        real_t *v = vectors + i * dim;

        for (uint_t j = 0; j < i; ++j) {

            real_t *u = vectors + j * dim;
            real_t r = dot(v, u, dim);

            for (uint_t k = 0; k < dim; ++k)
                v[k] = v[k] - r * u[k];

            R[IDX(j, i, dim)] = r;
        }
        R[IDX(i, i, dim)] = normalize(v, dim);
    }
}

int DLL_EXPORT extrapolate(real_t *vectors, uint_t dim, uint_t order, enum EXTR_METHOD method)
{
    real_t *const v0 = vectors;

    /* Compute differences */
    for (uint_t i = order - 1; i > 0; --i) {
        real_t *v = vectors + i * dim;
        real_t *u = v - dim;
        for (uint_t k = 0; k < dim; ++k)
            v[k] = v[k] - u[k];
    }
    vectors = vectors + dim;
    order = order - 1;

    /* QR decomposition */
    real_t R[order * order];
    mgs(vectors, R, dim, order);

    /* Compute coefficients */
    real_t gamma[order];
    switch (method) {
        case MPE:	mpe(R, gamma, order); break;
        case RRE:  	rre(R, gamma, order); break;
        default:    return -1;
    }

    /* xi = 1 - cumsum(gamma) */
    for (uint_t i = 1; i < order; ++i)
        gamma[i] = gamma[i] + gamma[i-1];
    real_t xi[order];
    for (uint_t i = 0; i < order; ++i)
        xi[i] = 1 - gamma[i];

    /* eta = R * xi */
    real_t eta[order];
    for (uint_t i = 0; i < order; ++i) {
        eta[i] = 0;
        for (uint_t k = 0; k < order; ++k)
            eta[i] = eta[i] + R[IDX(i, k, order)] * xi[k];
    }

    /* Compute the solution: v0 + Q * eta */
    for (uint_t i = 0; i < order; ++i) {
        real_t *q = vectors + i * dim;
        for (uint_t k = 0; k < dim; ++k)
            v0[k] = v0[k] + q[k] * eta[i];
    }
    return 0;
}
