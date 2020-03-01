#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>

int64_t cpucycles(void);
//½ÇÇà ÇÔ¼ö
void LtoR_Jacobian_Comb();
void LtoR_Jacobian_wNAF();
//
void addition(uint32_t* opA, uint32_t* opB, uint32_t* PF);
void addition_for_bitshift(uint32_t* opA, uint32_t* opB, uint32_t* PF);
void subtraction(uint32_t* opA, uint32_t* opB, uint32_t* PF);
void multiplication_os_64bit(uint32_t* A, uint32_t* B, uint32_t* C);
void square_64bit(uint32_t* A, uint32_t* C);
void Fast_Reduction_p256(uint32_t* c, uint32_t* result);
void binary_inversion(uint32_t* a);
//
void ECDBL_Jacobian(struct Jacobian* P1, struct Jacobian* P3);
int ECADD_Jacobian(struct Jacobian* P1, struct Affin* Q2, struct Jacobian* P3);
void change_affin(struct Jacobian* P);
//
void one_bit_shift(uint32_t* A, int i);
int equal_zero(uint32_t* A);
void array8_copy(uint32_t* A, uint32_t* B);
//
int not_equal_one_AB(uint32_t* A, uint32_t* B);
int greater_or_equal_AB(uint32_t* A, uint32_t* B);
void bit32_change_bit64(uint32_t* A, uint64_t* B);
void array16_copy(uint32_t* A, uint32_t* B);
int greater_or_equal_one(uint32_t* k);
//
void wNAF(uint32_t* k, char* k_i);
void wNAF_sm_version(char* k_i, struct Jacobian* Q);
void discriminate_mod16(char* k_i, int i, struct Jacobian* P);
//
void Fixed_Comb_Bit_table(char(*table)[32], uint32_t* k);
void Fixed_Comb_P_table(char(*table)[32], int j, struct Affin* P);
void ECLtoR_Fixed_Comb(char(*table)[32], struct Jacobian* Q);


//Affin ÁÂÇ¥°è
struct Affin {
	uint32_t X[16];
	uint32_t Y[16];
};

//Jacobian ÁÂÇ¥°è
struct Jacobian {
	uint32_t X[16];
	uint32_t Y[16];
	uint32_t Z[16];
};

