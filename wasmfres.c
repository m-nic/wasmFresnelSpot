#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3.h"
#include <complex.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#define ind_sq(i, j) i * N + j
#define ind_sq_sh(i, j) ((N/2 + i) % N) * N + ((N/2 + j) % N)

void printfftMat(fftw_complex *out, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("(%lf, %lf .j) ", out[ind_sq(i, j)][0], out[ind_sq(i, j)][1]);
        }
        printf("\n");
    }
}

void printMatJson(double *out, int N) {
    printf("[\n");
    for (int i = 0; i < N; i++) {
        printf("\t[");
        for (int j = 0; j < N; j++) {
            printf("%lf", out[ind_sq(i, j)]);
            if (j < N - 1) {
                printf(", ");
            }
        }

        printf("]");

        if (i < N - 1) {
            printf(",");
        }

        printf("\n");
    }
    printf("]");
}

fftw_complex *fft2d(fftw_complex *in, int N) {
    fftw_plan forward;

    fftw_complex *out;
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N * N);

    forward = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(forward);
    fftw_destroy_plan(forward);

    return out;
}

int hasRun = 0;
complex double *t2;
complex double *t3;
complex double *t4;
float lastZ = -1;

double *fresnel(int *M, int N, float z, float wv) {
    float k = 2 * M_PI / wv;

    if (hasRun == 0 || lastZ != z) {
        complex double t1 = cexp(I * k * z) / (I * wv * z);
        t2 = malloc(sizeof(complex double) * N * N);
        t3 = malloc(sizeof(complex double) * N * N);
        t4 = malloc(sizeof(complex double) * N * N);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int x0 = (-(N - 1) / 2) + i;
                int y0 = (-(N - 1) / 2) + j;
                int x1 = i + 1;
                int y1 = j + 1;

                t2[ind_sq(i, j)] = cexp(I * k * (pow(x0, 2) + pow(y0, 2)) / 2 * z);
                t3[ind_sq(i, j)] = cexp(I * k * (pow(x1, 2) + pow(y1, 2)) / 2 * z);

                t4[ind_sq(i, j)] = t1 * t2[ind_sq(i, j)];
            }
        }

        lastZ = z;
        hasRun = 1;
    }

    fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N * N);

    // dot mul M * t3 & shift
    complex double Mt3;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Mt3 = M[ind_sq(i, j)] * t3[ind_sq(i, j)];

            in[ind_sq_sh(i, j)][0] = creal(Mt3);
            in[ind_sq_sh(i, j)][1] = cimag(Mt3);
        }
    }
    //  FFT 2
    fftw_complex *fttRez = fft2d(in, N);
    double *FF = malloc(sizeof(double) * N * N);

    complex double U;
    double max = -INFINITY;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            U = t4[ind_sq(i, j)] * (fttRez[ind_sq_sh(i, j)][0] + fttRez[ind_sq_sh(i, j)][1] * I);
            FF[ind_sq(i, j)] = atan2(cimag(U), creal(U));

            if (max < FF[ind_sq(i, j)]) {
                max = FF[ind_sq(i, j)];
            }
        }
    }

    for (int i = 0; i < N * N; i++) {
        FF[i] = FF[i] / max;
    }

    // cleanup
//    free(t2);
//    free(t3);
//    free(t4);

    fftw_free(in);
    fftw_free(fttRez);

    return FF;
}

double *fresnelCircle(int xc, int yc, int r, float z, float wv, int N) {
    int *input = (int *) fftw_malloc(sizeof(int) * N * N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            input[ind_sq(i, j)] = sqrt(pow(i - xc, 2) + pow(j - yc, 2)) < r;
        }
    }

    double *out = fresnel(input, N, z, wv);

    fftw_free(input);

    return out;
}

int main() {
    int N = 512;

    int r = 50;
    int xc = N / 2;
    int yc = N / 2;

    double *out = fresnelCircle(xc, yc, r, 1e-9, 6.33e-6, N);

//    printMatJson(out, N);
    free(out);

    return 0;
}