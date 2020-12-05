#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "core.h"
#include "imcore.h"

// struct to keep four neighbours of the (i,j) and itself
struct sparse_element
{
    int index[5];
    float value[5];
};

// set the circular poission matrix elements
void fill_poisson_matrix(struct sparse_element *out, uint32_t r, uint32_t c, uint32_t nrows, uint32_t ncols, float beta)
{
    uint32_t r0 = (nrows + r - 1) % nrows;
    uint32_t r1 = (nrows + r + 0) % nrows;
    uint32_t r2 = (nrows + r + 1) % nrows;

    uint32_t c0 = (ncols + c - 1) % ncols;
    uint32_t c1 = (ncols + c + 0) % ncols;
    uint32_t c2 = (ncols + c + 1) % ncols;

    out->index[0] = ncols * r1 + c1;
    out->value[0] = 1 + 4 * beta;

    out->index[1] = ncols * r0 + c1;
    out->value[1] = -beta;

    out->index[2] = ncols * r2 + c1;
    out->value[2] = -beta;

    out->index[3] = ncols * r1 + c0;
    out->value[3] = -beta;

    out->index[4] = ncols * r1 + c2;
    out->value[4] = -beta;
}
                
// solve I = inv(P)*S
double sparse_solver(struct sparse_element *A, float *b, float *x, uint32_t length)
{
    float w = 1.5;
    uint32_t Max_Iter = 500;

    // fill x with zeros
    memset(x, length * sizeof(float), 0);

    // do SOR iterations
    double error = 0;
    uint32_t k, i, n;
    for (k = 0; k < Max_Iter; k++)
    {
        // do calculation for each pixel inside omega
        for (i = 0; i < length; i++)
        {
            // find the additional values coming from the border pixels
            float sum = 0;
            for(n = 1; n < 5; n++)
            {
                sum += A[i].value[n] * x[A[i].index[n]];
            }
            
            // do single SOR iteration
            x[i] = (1 - w) * x[i] + (w / A[i].value[0]) * (b[i] - sum);
        }

        // compute error
        error = 0;
        for (i = 0; i < length; i++)
        {
            // find the additional values coming from the border pixels
            float be = 0;
            for(n = 0; n < 5; n++)
            {
                be += A[i].value[n] * x[A[i].index[n]];
            }
            error += square(be - b[i]);
        }
        
        // find the mse error
        error = sqrt(error / length);

        // if the error is smaller than a threshold, finish the iterations
        if(error < 1e-5)
        {
            break;
        }
    }

    // return the error
    return error;
}

float gradient8y(matrix_t *input, uint32_t r, uint32_t c, uint32_t k)
{
    uint32_t r1 = (rows(input) + r + 0) % rows(input);
    uint32_t r2 = (rows(input) + r + 1) % rows(input);

    return (((float)atui8(input, r2, c, k)) - ((float)atui8(input, r1, c, k))) / 255.0f;
}

float gradient8x(matrix_t *input, uint32_t r, uint32_t c, uint32_t k)
{
    uint32_t c1 = (cols(input) + c + 0) % cols(input);
    uint32_t c2 = (cols(input) + c + 1) % cols(input);

    return (((float)atui8(input, r, c2, k)) - ((float)atui8(input, r, c1, k))) / 255.0f;
}

float gradient32y(matrix_t *input, uint32_t r, uint32_t c, uint32_t k)
{
    uint32_t r1 = (rows(input) + r + 0) % rows(input);
    uint32_t r2 = (rows(input) + r - 1) % rows(input);
    
    return atf(input, r2, c, k) - atf(input, r1, c, k);
}

float gradient32x(matrix_t *input, uint32_t r, uint32_t c, uint32_t k)
{
    uint32_t c1 = (cols(input) + c + 0) % cols(input);
    uint32_t c2 = (cols(input) + c - 1) % cols(input);

    return atf(input, r, c2, k) - atf(input, r, c1, k);
}


// convert float array to uint8_t image
void float2matrix(float **values, matrix_t *out)
{
    uint8_t *img = mdata(out, 0);

    // go over the data and set values
    uint32_t r,c,k, idx = 0;
    for(r = 0; r < rows(out); r++) 
    {
        for(c = 0; c < cols(out); c++) 
        {
            for(k = 0; k < channels(out); k++)
            {
                atui8(out, r, c, k) = clamp(values[k][idx] * 255, 0, 255);
            }
            // increase the index
            idx++;
        }
    }
}

// this functions finds the x,y gradient of the input and removes the smaller gradients less than the given threshold
void find_gradient(matrix_t *input, matrix_t *gx, matrix_t *gy, float threshold)
{
    // resize outputs
    matrix_resize(gx, rows(input), cols(input), channels(input));
    matrix_resize(gy, rows(input), cols(input), channels(input));

    // compute the gradient in the inside of the image
    uint32_t r, c, k;
    for(r = 0; r < rows(input); r++) 
    {
        for(c = 0; c < cols(input); c++) 
        {
            for(k = 0; k < channels(input); k++)
            {
                // get the x and y gradients
                float dx = gradient8x(input, r,c,k);
                float dy = gradient8y(input, r,c,k);

                // check that the gradient is large enough
                if((square(dx) + square(dy)) > threshold)
                {
                    atf(gx, r, c, k) = dx;
                    atf(gy, r, c, k) = dy;
                }
                else
                {
                    atf(gx, r, c, k) = 0.0f;
                    atf(gy, r, c, k) = 0.0f;
                }
            }
        }
    }
    //done
}

// this method implements the algorithm presented in Image Smoothing via L0 Gradient Minimization
void L0Minimize(matrix_t *input, float lambda, float beta0, float betaMax, float kappa, matrix_t *output)
{
    // resize the output matrix
    matrix_resize(output, rows(input), cols(input), channels(input));

    // initialization
    matrix_copy(input, output);
    float beta = beta0;

    // create variables defined in paper
    matrix_t *h = matrix_create(float);
    matrix_t *v = matrix_create(float);

    // create auxialary variables for iterations
    uint32_t r,c,k, idx;
    float **xvalues = malloc(sizeof(float*) * channels(input));
    float **bvalues = malloc(sizeof(float*) * channels(input));
    for(k = 0; k < channels(input); k++)
    {
        xvalues[k] = malloc(sizeof(float) * rows(input) * cols(input));
        bvalues[k] = malloc(sizeof(float) * rows(input) * cols(input));
    }

    struct sparse_element *Avalues = malloc(sizeof(struct sparse_element) * rows(input) * cols(input));

    // iteration
    uint32_t iter = 0;
    while(beta < betaMax)
    {
        // find dx and dy (gradient of the S) and modify such that hp,vp = 0 if (hp^2+vp^2) <= threshold
        find_gradient(output, h, v, lambda / beta);
        
        // solve S_{i+1} using h_i, v_i
        // S_{i+1} = inv(1 + beta(d^2/dx^2 + d^2/dy^2)) * (input + beta(dh_i/ dx + dv_i/dy))
        // S_{i+1} = inv(A) * b where A = 1 + beta(d^2/dx^2 + d^2/dy^2) and b = input + beta(dh_i/ dx + dv_i/dy)
        idx = 0;
        for(r = 0; r < rows(input); r++) 
        {
            for(c = 0; c < cols(input); c++) 
            {
                // create matrix A in sparse n vector form
                fill_poisson_matrix(&Avalues[idx], r, c, rows(input), cols(input), beta);

                for(k = 0; k < channels(input); k++)
                {
                    // create matrix b in column vector form
                    // get the x and y gradients of h,v
                    float dhx = gradient32x(h, r,c,k);
                    float dvy = gradient32y(v, r,c,k);

                    // compute b
                    bvalues[k][idx] = ((float)atui8(input, r, c, k)) / 255.0f + beta * (dhx + dvy);
                }
                // increase the index
                idx++;
            }
        }

        // solve Ax = b using successive over relaxation
        double err = 0;
        for(k = 0; k < channels(input); k++)
        {
            err += sparse_solver(Avalues, bvalues[k], xvalues[k], rows(input) * cols(input));
        }
        
        // print iteration number
        printf("Iteration[%03d]: %5.4f\r", iter++, err);
        fflush(stdout);

        // use xvalues to create S_{i+1}
        float2matrix(xvalues, output);

        // update beta
        beta = beta * kappa;
    }

    // deallocate unused memories
    for(k = 0; k < channels(input); k++)
    {
        free(xvalues[k]);
        free(bvalues[k]);
    }

    free(xvalues);
    free(bvalues);
    free(Avalues);

    matrix_free(&h);
    matrix_free(&v);
}

int main(int argc, unsigned char *argv[]) 
{
    if(argc != 2)
    {
        printf("call executable with a filename!\n");
        return -1;
    }

    matrix_t *img = imread(argv[1]);
    matrix_t *output = matrix_create(uint8_t);

    // set the L0 minimization paramaters
    float lambda = 0.03f;
    float kappa = 2.0f;
    float beta0 = 2.0f * lambda;
    float betaMax = 1e5f;

    // create L0 regularized image
    L0Minimize(img, lambda, beta0, betaMax, kappa, output);

    // write the results
    imwrite(output, "test_result.bmp");

    return 0;
}