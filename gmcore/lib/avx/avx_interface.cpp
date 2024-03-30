#include <immintrin.h>
// #include <iostream>
#include <stdio.h>

// void avx512_dot_product(double *a, double *b, int n, double *result);

extern "C" {
    void avx512_dot_product_(double *a, double *b, int *n, double *result);
}

void avx512_dot_product_(double *a, double *b, int *n, double *result) {
    int i;
    __m512d sum = _mm512_setzero_pd();
    // std::cout << n << std::endl;
    // printf("%d\n", n[0]);
    // printf("%d\n", a);
    for (i = 0; i < n[0]; i += 8) {
    // for (i = 0; i < 64; i += 8) {
        __m512d va = _mm512_loadu_pd(&a[i]);
        __m512d vb = _mm512_loadu_pd(&b[i]);
        sum = _mm512_add_pd(sum, _mm512_mul_pd(va, vb));
    }
    *result = _mm512_reduce_add_pd(sum);
}

// int main() {
//     return 0;
// }