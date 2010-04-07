#include <math.h>

typedef double real_t;
typedef unsigned long int size_t;

enum {MPE, RRE};

real_t dot(const real_t *v, const real_t *u, size_t dim)
{
    real_t res = 0;
    for (size_t k = 0; k < dim; ++k)
        res = res + v[k]*u[k];
    return res;
}

real_t normalize(real_t *v, size_t dim)
{
    real_t norm = sqrt(dot(v, v, dim));
    
    real_t inv_norm = 1.0 / norm;    
    for (size_t k = 0; k < dim; ++k)
        v[k] = v[k] * inv_norm;
        
    return norm;
}

void mgs(real_t *vectors, real_t *R, size_t dim, size_t num)
{
    for (size_t k = 0; k < num * num; ++k)
        R[k] = 0;
        
    for (size_t i = 0; i < num; ++i) {
        real_t *v = vectors + i * dim;

        for (size_t j = 0; j < i; ++j) {
            
            real_t *u = vectors + j * dim;
            real_t r = dot(v, u, dim);
            
            for (size_t k = 0; k < dim; ++k)
                v[k] = v[k] - r * u[k];

            R[j*dim + i] = r;
        }
        R[i*dim + i] = normalize(v, dim);
    }
}

int extrapolate(real_t *vectors, size_t dim, size_t num, int method)
{
    size_t i;
    real_t *const v0 = vectors;

    /* Compute differences */
    for (i = num - 1; i > 0; --i) {
        real_t *v = vectors + i * dim;
        real_t *u = v - dim;
        for (size_t k = 0; k < dim; ++k)
            v[k] = v[k] - u[k];
    }
    vectors = vectors + dim;
    num = num - 1;

    /* QR decomposition */
    real_t R[num * num];
    mgs(vectors, R, dim, num);
    
    /* Compute coefficients */
    switch (method) {
        case MPE:
            break;
        case RRE:
            break;
        default:
            return 0;
    }
    real_t eta[num];
    
    /* Compute the solution */
    for (size_t k = 0; k < num; ++k) {
        real_t *u = vectors + k * dim;
        int i;
        for (i = 0; i < dim; ++i) {        
            v0[i] = v0[i] + eta[k] * u[i];
        }
    }
    return 1;
}
