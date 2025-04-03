#include "../include/matrix.h"
#include "../include/param.h"
#include <string.h>

void scloudplus_add(uint16_t *in0, uint16_t *in1, int len, uint16_t *out) {
    for (int i = 0; i < len; i++) {
        out[i] = (in0[i] + in1[i]) & 0xFFF;
    }
}

void scloudplus_sub(uint16_t *in0, uint16_t *in1, int len, uint16_t *out) {
    for (int i = 0; i < len; i++) {
        out[i] = (in0[i] - in1[i]) & 0xFFF;
    }
}

void scloudplus_mul_cs(uint16_t *C, uint16_t *S,
                       uint16_t *out) {

    memset(out, 0, scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));

    for (int i = 0; i < scloudplus_mbar; i++) {
        uint16_t *out_row = &out[i * scloudplus_nbar];

        for (int k = 0; k < scloudplus_n; k++) {
            const uint16_t c_val = C[i * scloudplus_n + k];

            for (int j = 0; j < scloudplus_nbar; j++) {
                const uint16_t s_val = S[j * scloudplus_n + k];
                out_row[j] += c_val * s_val;
            }
        }
    }
}

void scloudplus_mul_add_sb_e(const uint16_t *S,
                             const uint16_t *B,
                             const uint16_t *E,
                             uint16_t *out) {
    memcpy(out, E, scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));

    for (int i = 0; i < scloudplus_mbar; i++) {
        const uint16_t *S_row = &S[i * scloudplus_m];

        for (int k = 0; k < scloudplus_m; k++) {
            const uint16_t s_val = S_row[k];
            const uint16_t *B_row = &B[k * scloudplus_nbar];
            uint16_t *out_row = &out[i * scloudplus_nbar];

            int j = 0;
            for (; j < (scloudplus_nbar & ~3); j += 4) {
                out_row[j] += s_val * B_row[j];
                out_row[j + 1] += s_val * B_row[j + 1];
                out_row[j + 2] += s_val * B_row[j + 2];
                out_row[j + 3] += s_val * B_row[j + 3];
            }

            for (; j < scloudplus_nbar; j++) {
                out_row[j] += s_val * B_row[j];
            }
        }
    }
}
