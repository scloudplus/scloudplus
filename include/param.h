#ifndef _SCLOUDPLUS_PARAM_H_
#define _SCLOUDPLUS_PARAM_H_

#ifndef scloudplus_l
#define scloudplus_l 128
#endif
#if (scloudplus_l == 128)
#define scloudplus_ss 16
#define scloudplus_m 600
#define scloudplus_n 600
#define scloudplus_mbar 8
#define scloudplus_nbar 8
#define scloudplus_logq 12
#define scloudplus_logq1 9
#define scloudplus_logq2 7
#define scloudplus_h1 150
#define scloudplus_h2 150
#define scloudplus_eta1 7
#define scloudplus_eta2 7
#define scloudplus_mu 64
#define scloudplus_tau 3
#define scloudplus_subm 2
#define scloudplus_block_number 75
#define scloudplus_block_size 4
#define scloudplus_block_rowlen 300
#define scloudplus_c1 5400 // MBAR*N*LOGQ1/sizeof(uint8_t)
#define scloudplus_c2 56
#define scloudplus_ctx 5456
#define scloudplus_pk 7216
#define scloudplus_pke_sk 1200
#define scloudplus_kem_sk 8480
#define scloudplus_m2 360000
#define scloudplus_m3 216000000
#define scloudplus_n2 360000
#define scloudplus_n3 216000000
#define scloudplus_mnin 679
#define scloudplus_mnout 582
#elif (scloudplus_l == 192)
#define scloudplus_ss 24
#define scloudplus_m 928
#define scloudplus_n 896
#define scloudplus_mbar 8
#define scloudplus_nbar 8
#define scloudplus_logq 12
#define scloudplus_logq1 12
#define scloudplus_logq2 10
#define scloudplus_h1 224
#define scloudplus_h2 232
#define scloudplus_eta1 2
#define scloudplus_eta2 1
#define scloudplus_mu 96
#define scloudplus_tau 4
#define scloudplus_subm 2
#define scloudplus_block_number 112
#define scloudplus_block_size 4
#define scloudplus_block_rowlen 448
#define scloudplus_c1 10752
#define scloudplus_c2 80
#define scloudplus_ctx 10832
#define scloudplus_pk 11152
#define scloudplus_pke_sk 1792
#define scloudplus_kem_sk 13008
#define scloudplus_mnin 671
#define scloudplus_mnout 488
#elif (scloudplus_l == 256)
#define scloudplus_ss 32
#define scloudplus_m 1136
#define scloudplus_n 1120
#define scloudplus_mbar 12
#define scloudplus_nbar 11
#define scloudplus_logq 12
#define scloudplus_logq1 10
#define scloudplus_logq2 7
#define scloudplus_h1 280
#define scloudplus_h2 284
#define scloudplus_eta1 3
#define scloudplus_eta2 2
#define scloudplus_mu 64
#define scloudplus_tau 3
#define scloudplus_subm 4
#define scloudplus_block_number 140
#define scloudplus_block_size 4
#define scloudplus_block_rowlen 560
#define scloudplus_c1 16800
#define scloudplus_c2 116
#define scloudplus_ctx 16916
#define scloudplus_pk 18760
#define scloudplus_pke_sk 3080
#define scloudplus_kem_sk 21904
#define scloudplus_n2 1254400
#define scloudplus_n3 1404928000
#define scloudplus_n4 1573519360000
#define scloudplus_n5 1762341683200000
#define scloudplus_m2 1290496
#define scloudplus_m3 1466003456
#define scloudplus_m4 1665379926016
#define scloudplus_m5 1891871595954176
#define scloudplus_mnin 680
#define scloudplus_mnout 530
#else
#error \
	"Please define one of scloudplus_l_128, scloudplus_l_192, or scloudplus_l_256"
#endif
#endif // _SCLOUDPLUS_PARAM_H_
