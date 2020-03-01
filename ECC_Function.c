#include "ECC.h"

int64_t cpucycles(void)
{
	return __rdtsc();
}

//main함수에 이 함수만 추가시면 ECLtoR_Fixed_Comb_version을 실행시킬수 있다.
void LtoR_Jacobian_Comb()
{
	struct Jacobian Q;
	struct Jacobian* Pointer_Q = &Q;

	uint32_t k[8] = { 0, };
	char table[8][32] = { 0, };
	int count = 0;
	int i = 0;
	unsigned long long cycles = 0, cycles1 = 0, cycles2 = 0, totalcycle = 0;


	FILE* file1 = fopen("TV_Scalar.txt", "r");
	FILE* file3 = fopen("TV_SM.txt", "w");
	while (count < 10000)
	{
		count++;
		printf("%d \n", count);

		//input
		for (i = 0; i < 8; i++)
		{
			fscanf(file1, "%08X", &k[7 - i]);
		}

		//cycle
		cycles1 = cpucycles();

		//중간 과정
		Fixed_Comb_Bit_table(table, k);
		ECLtoR_Fixed_Comb(table, Pointer_Q);

		//cycle
		cycles2 = cpucycles();
		cycles = cycles2 - cycles1;
		totalcycle += cycles;

		//output
		for (i = 0; i < 8; i++)
			fprintf(file3, "%08X", Q.X[7 - i]);
		fprintf(file3, "\n");
		for (i = 0; i < 8; i++)
			fprintf(file3, "%08X", Q.Y[7 - i]);
		fprintf(file3, "\n\n");
	}
	fclose(file1);
	fclose(file3);

	printf(" totalcycles/10000 = %10lld", totalcycle / 10000);

}

//main함수에 이 함수만 추가시면 wNAF_sm_version을 실행시킬수 있다.
void LtoR_Jacobian_wNAF()
{
	unsigned long long cycles = 0, cycles1 = 0, cycles2 = 0, totalcycle = 0;

	struct Jacobian Q;
	struct Jacobian* Pointer_Q = &Q;

	uint32_t k[8] = { 0, };
	char k_i[257] = { 0, };
	int count = 0;
	int i = 0;


	FILE* file1 = fopen("TV_Scalar.txt", "r");
	FILE* file3 = fopen("TV_SM.txt", "w");
	while (count < 10000)
	{
		count++;
		printf("%d \n", count);
		
		//input
		for (i = 0; i < 8; i++)
		{
			fscanf(file1, "%08X", &k[7 - i]);
		}
		//cycle
		cycles1 = cpucycles();
		
		//중간 과정
		wNAF(k, k_i);
		wNAF_sm_version(k_i, &Q);
		
		//cycle
		cycles2 = cpucycles();
		cycles = cycles2 - cycles1;
		totalcycle += cycles;

		//output
		for (i = 0; i < 8; i++)
			fprintf(file3, "%08X", Q.X[7 - i]);
		fprintf(file3, "\n");
		for (i = 0; i < 8; i++)
			fprintf(file3, "%08X", Q.Y[7 - i]);
		fprintf(file3, "\n\n");
	}
	fclose(file1);
	fclose(file3);

	printf(" totalcycles/10000 = %10lld", totalcycle / 10000);
}

//input : opA, opB  output : opA + opB = PF[8]
void addition(uint32_t* opA, uint32_t* opB, uint32_t* PF)
{
	uint32_t carry[9] = { 0, };
	uint32_t borrow[9] = { 0, };
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };
	uint32_t xr[1] = { 0xFFFFFFFF };
	uint32_t one[9] = { 0x00000001,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000 };
	int i = 0;

	//result 초기화
	for (i = 0; i < 9; i++)
	{
		PF[i] = 0;
	}
		
	for (i = 0; i < 8; i++)
	{
		PF[i] = opA[i] + opB[i];

		if (PF[i] < opA[i])
		{
			carry[i + 1] = 1;
		}

		PF[i] += carry[i];

		if (PF[i] < carry[i])
		{
			carry[i + 1] += 1;
		}

		if (carry[8] == 1)
		{
			PF[8] = 1;
		}
	}

	if (carry[8] == 1)
	{

		for (i = 0; i < 8; i++)//borrow 1step
		{
			if (PF[i] < P256[i])
			{
				borrow[i + 1] = 1;
			}
		}
		for (i = 0; i < 8; i++)//borrow 2step
		{
			if (PF[i] == P256[i] && borrow[i] == 1)
			{
				borrow[i + 1] = 1;
			}
		}
		for (i = 0; i < 8; i++)
		{
			if (borrow[i + 1] == 1)
			{
				PF[i] = PF[i] + 0xffffffff - P256[i];
				PF[i] = PF[i] + one[0];
				PF[i] = PF[i] - borrow[i];
			}
			else
			{
				PF[i] = PF[i] - P256[i] - borrow[i];
			}
		}
	}
	//result : PF[7~0]
}

void addition_for_bitshift(uint32_t *opA, uint32_t* opB, uint32_t* PF)
{
	int i;
	uint32_t carry[9] = { 0x00000000, };
	uint32_t borrow[9] = { 0x00000000, };
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };//{0x00000000 ,0xFFFFFFFF, 0x00000001,0x00000000, 0x00000000, 0x00000000,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF}
	uint32_t xr[1] = { 0xFFFFFFFF };
	uint32_t one[9] = { 0x00000001,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000 };

	for (i = 0; i < 9; i++)
	{
		PF[i] = 0;
		
	}

	for (i = 0; i < 8; i++)
	{
		PF[i] = opA[i] + opB[i];

		if (PF[i] < opA[i])
		{
			carry[i + 1] = 1;
		}

		PF[i] += carry[i];

		if (PF[i] < carry[i])
		{
			carry[i + 1] += 1;
		}

		if (carry[8] == 1)
		{
			PF[8] = 1;
		}
	}
	//result : PF[8~0] (PF[8]에 CARRY 1 이 발생되있을수 있음

}

//input : opA, opB  output : opA - opB = PF[8]
void subtraction(uint32_t* opA, uint32_t* opB, uint32_t* PF)
{
	int i;
	uint32_t carry[9] = { 0x00000000, };
	uint32_t borrow[9] = { 0x00000000, };
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };//{0x00000000 ,0xFFFFFFFF, 0x00000001,0x00000000, 0x00000000, 0x00000000,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF}
	uint32_t xr[1] = { 0xFFFFFFFF };
	uint32_t one[9] = { 0x00000001,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000 };
	
	//result 초기화
	for (i = 0; i < 9; i++)
	{
		PF[i] = 0;
	}

	for (i = 0; i < 8; i++)//1step
	{
		if (opA[i] < opB[i])
		{
			borrow[i + 1] = 1;
		}
	}

	for (i = 0; i < 8; i++)//2step
	{
		if (opA[i] == opB[i] && borrow[i] == 1)
		{
			borrow[i + 1] = 1;
		}
	}

	for (i = 0; i < 8; i++)
	{
		if (borrow[i + 1] == 1)
		{
			PF[i] = 0xffffffff - opB[i];
			PF[i] += opA[i];
			PF[i] += one[0];
			PF[i] -= borrow[i];
		}
		else
		{
			PF[i] = opA[i] - opB[i] - borrow[i];
		}
	}

	//버로우가 마지막에 1이면 -되는거니깐 P를 한번 더해줘야해 
	if (borrow[8] == 1)
	{
		for (i = 0; i < 8; i++)
		{
			PF[i] = PF[i] + P256[i];

			if (PF[i] < P256[i])
			{
				carry[i + 1] = 1;
			}

			PF[i] += carry[i];

			if (PF[i] < carry[i])
			{
				carry[i + 1] += 1;
			}

			if (carry[8] == 1)
			{
				PF[8] = 1;
			}
		}

	}
	//result : PF[7~0]
}

//input : opA, opB  output : opA*opB = C[16]
void multiplication_os_64bit(uint32_t* A, uint32_t* B, uint32_t* C)
{
	int i=0;
	int j=0;
	uint64_t UV=0;
	uint64_t opA[9] = { 0, };
	uint64_t opB[9] = { 0, };
	uint32_t U=0;
	uint32_t V=0;

	//32bit->64bit
	bit32_change_bit64(A, opA);
	bit32_change_bit64(B, opB);

	for (i = 0; i < 15; i++)
	{
		C[i] = 0;
	}

	for (i = 0; i < 8; i++)
	{
		U = 0;

		for (j = 0; j < 8; j++)
		{
			
			UV = C[i + j] + opA[j] * opB[i] + U;

			U = ((UV >> 32) & 0xffffffff);
			V = (UV & 0xffffffff);

			C[i + j] = V;
		}
		C[i + 8] = U;
	}
	//result : C[15~0]
}

//input : A  output : opA^2 = C[16]
void square_64bit(uint32_t* A, uint32_t* C)
{
	uint64_t UV = 0, UVCP = 0;
	uint64_t opA[16] = { 0, };
	uint32_t R0 = 0, R1 = 0, R2 = 0;
	uint32_t U = 0;
	uint32_t V = 0;

	int i, j, k;
	int carry = 0;

	bit32_change_bit64(A, opA);

	for (i = 0; i < 15; i++)
	{
		C[i] = 0;
	}

	for (k = 0; k < 15; k++)
	{
		for (i = k, j = 0; i >= j; i--, j++)
		{
			while (i >= 8)
			{
				i--;
				j++;
			}
			if (i == j)
			{
				UV = opA[i] * opA[j];
				//V추출
				V = (UV & 0xffffffff);
				//carry값 발생하지 32bit=32bit+32bit 
				R0 += V;
				if (R0 < V)
				{
					carry = 1;
				}
				//U추출하고 
				U = ((UV >> 32) & 0xffffffff);

				R1 = R1 + U + carry;

				//carry 초기화 
				carry = 0;

				if (R1 < U)
				{
					carry = 1;
				}
				R2 = R2 + carry;
				//carry 초기화2 
				carry = 0;

				break;
			}

			UV = opA[i] * opA[j];
			UVCP = UV;
			UV = UV + UVCP;
			if (UV < UVCP)
			{
				carry = 1;
			}
			R2 = R2 + carry;
			//carry 초기화2 
			carry = 0;

			//V추출하고 
			V = (UV & 0xffffffff);
			//carry값 발생하지 32bit=32bit+32bit 
			R0 = R0 + V;
			if (R0 < V)
			{
				carry = 1;
			}
			//U추출하고 
			U = ((UV >> 32) & 0xffffffff);

			R1 = R1 + U + carry;
			//carry 초기화 
			carry = 0;

			if (R1 < U)
			{
				carry = 1;
			}
			R2 = R2 + carry;
			//carry 초기화2 
			carry = 0;
		}
		C[k] = R0; R0 = R1; R1 = R2; R2 = 0;
	}
	C[15] = R0; //바로 윗줄 R0=R1때문에 
	
	//result : C[15~0]
}

//input : C[15~0]  output : (reduction p256)temp2[7~0]
void Fast_Reduction_p256(uint32_t *c, uint32_t* result)
{
	int i = 0;
	int j = 0;
	uint32_t s1[8], s2[8], s3[8], s4[8], s5[8], s6[8], s7[8], s8[8], s9[8];
	uint32_t temp1[9] = { 0, };
	
	for (i = 0; i < 16; i++)
	{
		result[i] = 0;
	}

	//s1
	for (i = 0; i < 8; i++)
	{
		s1[i] = c[i];
	}

	//s2
	s2[7] = c[15]; s2[6] = c[14]; s2[5] = c[13]; s2[4] = c[12]; s2[3] = c[11]; s2[2] = 0; s2[1] = 0; s2[0] = 0;

	//s3
	s3[7] = 0; s3[6] = c[15]; s3[5] = c[14]; s3[4] = c[13]; s3[3] = c[12]; s3[2] = 0; s3[1] = 0; s3[0] = 0;

	//s4
	s4[7] = c[15]; s4[6] = c[14]; s4[5] = 0; s4[4] = 0; s4[3] = 0; s4[2] = c[10]; s4[1] = c[9]; s4[0] = c[8];

	//s5
	s5[7] = c[8]; s5[6] = c[13]; s5[5] = c[15]; s5[4] = c[14]; s5[3] = c[13]; s5[2] = c[11]; s5[1] = c[10]; s5[0] = c[9];

	//s6
	s6[7] = c[10]; s6[6] = c[8]; s6[5] = 0; s6[4] = 0; s6[3] = 0; s6[2] = c[13]; s6[1] = c[12]; s6[0] = c[11];

	//s7
	s7[7] = c[11]; s7[6] = c[9]; s7[5] = 0; s7[4] = 0; s7[3] = c[15]; s7[2] = c[14]; s7[1] = c[13]; s7[0] = c[12];

	//s8
	s8[7] = c[12]; s8[6] = 0; s8[5] = c[10]; s8[4] = c[9]; s8[3] = c[8]; s8[2] = c[15]; s8[1] = c[14]; s8[0] = c[13];

	//s9
	s9[7] = c[13]; s9[6] = 0; s9[5] = c[11]; s9[4] = c[10]; s9[3] = c[9]; s9[2] = 0; s9[1] = c[15]; s9[0] = c[14];

	//s1+2(s2)
	addition(s1, s2, temp1);
	addition(temp1, s2, result);
	//s1+2(s3)
	addition(result, s3, temp1);
	addition(temp1, s3, result);
	//s1+s4
	addition(result, s4, temp1);
	//s1+s5
	addition(temp1, s5, result);
	//s1-s6
	subtraction(result, s6, temp1);
	//s1-s7
	subtraction(temp1, s7, result);
	//s1-s8
	subtraction(result, s8, temp1);
	//s1-s9
	subtraction(temp1, s9, result);

	//result : temp2[7~0]
}

//input : (A[8]and i=7) or (A[9] and i=8) , output : one bit(>>1) shift
void one_bit_shift(uint32_t* A, int i)
{
	uint32_t bit_shift[9] = { 0, };
	int j = 0;
	
	for (j=i ; j > 0; j--)
	{
		bit_shift[j] = A[j] & 0x01;
		bit_shift[j] = bit_shift[j] << 31;
	}
	for (j = i; j >= 0; j--)
	{
		A[j] = A[j] >> 1;
	}
	
	for (j = i; j > 0; j--)
	{
		A[j - 1] = A[j - 1] ^ bit_shift[j];
	}
	
	
}

//input : P=(X1,Y1,Z1)  output : 2P=(X3,Y3,Z3)
void ECDBL_Jacobian(struct Jacobian *P1, struct Jacobian* P3)
{
	int i = 0;
	uint32_t temp[16] = { 0, };
	uint32_t temp2[16] = { 0, };
	uint32_t T1[16] = { 0, };
	uint32_t T2[16] = { 0, };
	uint32_t T3[16] = { 0, };
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };
	
	//결과값 0으로 초기화
	uint32_t ZERO[16] = { 0, };
	array16_copy(ZERO, P3->X);
	array16_copy(ZERO, P3->Y);
	array16_copy(ZERO, P3->Z);

	//2step T1<-Z1^2
	square_64bit(P1->Z, temp);
	Fast_Reduction_p256(temp, T1);
	
	//3step
	subtraction(P1->X, T1, T2);
	
	//4step
	array8_copy(T1, temp);
	addition(P1->X, temp ,T1);

	//5step
	array8_copy(T2, temp);
	multiplication_os_64bit(temp, T1, T2);
	Fast_Reduction_p256(T2, temp);
	array8_copy(temp, T2);
	
	//6step
	addition(T2, T2, temp);
	addition(temp, T2, temp2);
	array8_copy(temp2, T2);

	//7step
	addition(P1->Y, P1->Y, P3->Y);

	//8step
	multiplication_os_64bit(P3->Y, P1->Z, P3->Z);
	Fast_Reduction_p256(P3->Z, temp);
	array8_copy(temp, P3->Z);

	//9step
	square_64bit(P3->Y, temp);
	Fast_Reduction_p256(temp, P3->Y);

	//10step
	multiplication_os_64bit(P3->Y, P1->X, temp);
	Fast_Reduction_p256(temp, T3);
	
	//11step
	square_64bit(P3->Y, temp);
	Fast_Reduction_p256(temp, P3->Y);
	
	//12step 중요
	if (((P3->Y[0]) & 0x01) == 0)
	{
		one_bit_shift(P3->Y, 7);
	}
	//else y3<-(y3+p)/2
	else
	{
		//y3+P256 (감산 안해줘)
		addition_for_bitshift(P3->Y, P256, temp);
		one_bit_shift(temp, 8);
		array8_copy(temp, P3->Y);
	}	

	//13step
	square_64bit(T2, temp);
	Fast_Reduction_p256(temp, P3->X);

	//14step
	addition(T3, T3, T1);

	//15step
	array8_copy(P3->X, temp);
	subtraction(temp, T1, P3->X);
		
	//16step
	subtraction(T3, P3->X, T1);

	//17step
	array8_copy(T1, temp);
	multiplication_os_64bit(temp, T2, T1);
	Fast_Reduction_p256(T1, temp);
	array8_copy(temp, T1);

	//18step
	array8_copy(P3->Y, temp);
	subtraction(T1, temp, P3->Y);

	// result : 2P = (x3,y3,z3)

}

//input : P1 is Jacobian, Q2 is Affin  , output : P3 is Jacobian
int ECADD_Jacobian(struct Jacobian* P1, struct Affin* Q2, struct Jacobian* P3)
{
	int i=0;
	uint32_t temp[16] = { 0, };
	uint32_t temp2[16] = { 0, };
	uint32_t T1[16] = { 0, };
	uint32_t T2[16] = { 0, };
	uint32_t T3[16] = { 0, };
	uint32_t T4[16] = { 0, };
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };

	//결과값 0으로 초기화
	uint32_t ZERO[16] = { 0, };
	array16_copy(ZERO, P3->X);
	array16_copy(ZERO, P3->Y);
	array16_copy(ZERO, P3->Z);

	//3STEP
	square_64bit(P1->Z, temp);
	Fast_Reduction_p256(temp, T1);

	//4STEP
	multiplication_os_64bit(T1, P1->Z, temp);
	Fast_Reduction_p256(temp, T2);

	//5STEP
	multiplication_os_64bit(T1, Q2->X, temp);
	Fast_Reduction_p256(temp, T1);

	//6step
	multiplication_os_64bit(T2, Q2->Y, temp);
	Fast_Reduction_p256(temp, T2);

	//7step
	array8_copy(T1, temp);
	subtraction(temp, P1->X, T1);

	//8step
	array8_copy(T2, temp);
	subtraction(temp, P1->Y, T2);

	//9step
	if ((equal_zero(T1)) == 0)
	{
		if ((equal_zero(T2)) == 0)
		{
			array8_copy(Q2->X, P1->X);
			array8_copy(Q2->Y, P1->Y);
			for (i = 8; i >0 ; i--)
			{
				P1->Z[i] = 0;
			}
			P1->Z[0] = 1;
			ECDBL_Jacobian(P1, P3); //T1 T2 BEFORE STEP9 => BOTH ZERO
			return 0;
		}
	}

	//10step
	multiplication_os_64bit(P1->Z, T1, temp);
	Fast_Reduction_p256(temp, P3->Z);

	//11step
	square_64bit(T1, temp);
	Fast_Reduction_p256(temp, T3);

	//12step
	multiplication_os_64bit(T3, T1, temp);
	Fast_Reduction_p256(temp, T4);

	//13step
	array8_copy(T3, temp2);
	multiplication_os_64bit(temp2, P1->X, temp);
	Fast_Reduction_p256(temp, T3);

	//14step
	addition(T3, T3, T1);

	//15step
	square_64bit(T2, temp);
	Fast_Reduction_p256(temp, P3->X);

	//16step
	array8_copy(P3->X, temp);
	subtraction(temp, T1, P3->X);
	
	//17step
	array8_copy(P3->X, temp);
	subtraction(temp, T4, P3->X);

	//18step
	array8_copy(T3, temp);
	subtraction(temp, P3->X, T3);

	//19step
	multiplication_os_64bit(T3, T2, temp);
	Fast_Reduction_p256(temp, T3);

	//20step
	multiplication_os_64bit(T4, P1->Y, temp);
	Fast_Reduction_p256(temp, T4);

	//21step
	subtraction(T3, T4, P3->Y);

	return 0;
	//return output : P3 = P+Q = (x3:y3:z3) in jacobian
}

void binary_inversion(uint32_t *a)
{
	uint32_t P256[9] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF, 0x00000000 };
	uint32_t u[9] = { 0, };
	uint32_t v[9] = { 0, };
	uint32_t x1[9] = { 0, };
	uint32_t x2[9] = { 0, };
	uint32_t temp[9] = { 0, };
	
	int i = 0;
	//1 Line
	array8_copy(a, u);
	array8_copy(P256, v);

	//2Line
	x1[0] = 1;
	x2[0] = 0;

	//3Line
	while (not_equal_one_AB(u, v) == -1)
	{
		//3.1 Line
		while ((u[0] & 0x01) == 0)
		{
			//u<-u/2
			one_bit_shift(u, 7);
			
			//If x1 is even, then x1<-x1/2
			if ((x1[0] & 0x01) == 0)
			{
				one_bit_shift(x1, 7);
			}
			//else x1<-(x1+p)/2
			else
			{
				//x1+P256(감산x)
				array8_copy(x1, temp);
				addition_for_bitshift(temp, P256, x1);
				
				//x1/2
				one_bit_shift(x1, 8);
			}
		}
		
		//3.2 Line
		while ((v[0] & 0x01) == 0)
		{
			//v<-v/2.
			one_bit_shift(v, 7);

			//if x2 is even, then x2<-x2/2
			if ((x2[0] & 0x01) == 0)
			{
				one_bit_shift(x2, 7);
			}
			//else x2<-(x2+p)/2
			else
			{
				//x2+P256
				array8_copy(x2, temp);
				addition_for_bitshift(temp, P256, x2);

				//x2/2
				one_bit_shift(x2, 8);
			}
		}

		//3.3Line
		if (greater_or_equal_AB(u, v) == 1)
		{
			//u<-u-v
			array8_copy(u, temp);
			subtraction(temp, v, u);
			//x1<-x1-x2
			array8_copy(x1, temp);
			subtraction(temp, x2, x1);
		}
		else
		{
			//v<-v-u
			array8_copy(v, temp);
			subtraction(temp, u, v);

			//x2<-x2-x1
			array8_copy(x2, temp);
			subtraction(temp, x1, x2);
		}
	}

	//4Line
	if (u[0] == 1)
	{
		array8_copy(x1, a);
	}
	else
	{
		array8_copy(x2, a);	
	}
	//result : invers of a
}

//A[8]배열이 0이면 0을 return , 1이면 1을 return
int equal_zero(uint32_t* A)
{
	int i=0;
	for (i = 7; i > 0; i--)
	{
		if (A[i] !=  0)
		{
			return 1;
		}
	}
	if (A[0] == 0)
	{
		return 0;
	}
	return 1;
}

//A[8]를 B[8]에 복붙한다.
void array8_copy(uint32_t* A, uint32_t* B)
{
	int i=0;
	//대상 초기화
	for (i = 0; i < 8; i++)
	{
		B[i] = 0;
	}

	for (i = 0; i < 8; i++)
	{
		B[i] = A[i];
	}
}
void array16_copy(uint32_t* A, uint32_t* B)
{
	int i = 0;
	//대상 초기화
	for (i = 0; i < 16; i++)
	{
		B[i] = 0;
	}

	for (i = 0; i < 16; i++)
	{
		B[i] = A[i];
	}
}

void change_affin(struct Jacobian* P)
{
	int i = 0;
	uint32_t temp[16] = { 0, };
	uint32_t temp2[16] = { 0, };

	array8_copy(P->Z, temp2);
	//Z^2
	square_64bit(P->Z, temp);
	Fast_Reduction_p256(temp, P->Z);
	//Z^3
	multiplication_os_64bit(P->Z, temp2, temp);
	Fast_Reduction_p256(temp, P->Z);
	//1/Z^3
	binary_inversion(P->Z);

	//Affin Y축
	multiplication_os_64bit(P->Y, P->Z, temp);
	Fast_Reduction_p256(temp, P->Y);

	multiplication_os_64bit(temp2, P->Z, temp);
	Fast_Reduction_p256(temp, P->Z);
	//Affin X축
	multiplication_os_64bit(P->X, P->Z, temp);
	Fast_Reduction_p256(temp, P->X);


}

//for binary_inversion : u!=1 and v!=1 일 경우 -1 return  아닐경우 1 return
int not_equal_one_AB(uint32_t *A, uint32_t *B)
{
	int i = 0;
	for (i = 7; i > 0; i--)
	{
		if ((A[i] != 0) && (B[i] != 0))
		{
			return -1;
		}

	}
	if ((A[0] != 1) && B[0] != 1)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

int greater_or_equal_AB(uint32_t* A, uint32_t* B)
{
	int i = 0;
	for (i = 7; i >= 0; i--)
	{
		if (A[i] > B[i])
		{
			return 1;
		}

		if (A[i] < B[i])
		{
			return -1;
		}
	}
	return 1;
}

void bit32_change_bit64(uint32_t *A, uint64_t *B)
{
	int i = 0;
	for (i = 0; i < 8; i++)
	{
		B[i] = A[i];
	}

}

void wNAF(uint32_t *k, char* k_i)
{
	//결과값 초기화
	int j = 0;
	for(j = 0; j < 257; j++)
	{
		k_i[j] = 0;
	}

	//1Line
	int i = 0;
	//2Line
	while (greater_or_equal_one(k) == 1)
	{
		//2.1Line
		if ((k[0] & 1) == 1)
		{
			k_i[i] = k[0] & 0xf;
			if (k_i[i] > 8)
			{
				k_i[i] -= 16;
			}
			k[0] = k[0] - k_i[i];
		}
		//2.2Line
		else
		{
			k_i[i] = 0;
		}
		
		//2.3Line
		one_bit_shift(k, 7);
		
		i++;
		
	}
}

void wNAF_sm_version(char *k_i, struct Jacobian* Q)
{
	struct Jacobian P;
	struct Jacobian* Pointer_P = &P;

	struct Jacobian temp_P, temp_Q;
	struct Jacobian* Pointer_tP = &temp_P, * Pointer_tQ = &temp_Q;

	struct Affin Q_Affin;
	struct Affin* Pointer_Q_Affin = &Q_Affin;
	
	int i = 0;
	int j = 0;
	int count = 0;

	//4Line
	for (i = 256; i >= 0; i--)
	{
		//4.1Line : Q<-2Q (ECDBL_Jacobian)
		if (count == 1)
		{
			array16_copy(Q->X, Pointer_tQ->X);
			array16_copy(Q->Y, Pointer_tQ->Y);
			array16_copy(Q->Z, Pointer_tQ->Z);

			ECDBL_Jacobian(Pointer_tQ, Q);
		}
		//4.2Line
		if (k_i[i] != 0)
		{
			//처음에는 ECADD_Jacobian 해줄 필요 없으니깐
			if (count == 0)
			{
				//iP 결정
				discriminate_mod16(k_i, i, Pointer_P);
				
				array16_copy(Pointer_P->X, Q->X);
				array16_copy(Pointer_P->Y, Q->Y);
				array16_copy(Pointer_P->Z, Q->Z);
				
				count=1;
			}
			//두번째부터는 ECADD_Jacobian 해주어야함
			else
			{
				//iP 결정
				discriminate_mod16(k_i, i, Pointer_P);
				//Q(Jacobian)을 Affin_Q로 change
				change_affin(Q);
				//Q를 struct Affin Q_Affin으로 대입
				array16_copy(Q->X, Pointer_Q_Affin->X);
				array16_copy(Q->Y, Pointer_Q_Affin->Y);
				
				ECADD_Jacobian(&P, Pointer_Q_Affin, Q);
			}
		}
		
	
	
	}
	//Q를 Affin으로 바꾸어주어야한다.
	change_affin(Q);
}


//판별 and table
void discriminate_mod16(char* k_i, int i, struct Jacobian *P)
{
	//struct Jacobian One_P, Three_P, Five_P, Seven_P, mOne_P, mThree_P, mFive_P, mSeven_P;
	//struct Jacobian* Pointer_One_P= &One_P, * Pointer_Three_P = &Three_P, * Pointer_Five_P = &Five_P, * Pointer_Seven_P = &Seven_P, * Pointer_mOne_P = &mOne_P, * Pointer_mThree_P = &mThree_P, * Pointer_mFive_P = &mFive_P, * Pointer_mSeven_P = &mSeven_P;
	
	//wNAF->table
	uint32_t One_P_x[16] = { 0xD898C296, 0xF4A13945, 0x2DEB33A0, 0x77037D81, 0x63A440F2, 0xF8BCE6E5, 0xE12C4247, 0x6B17D1F2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t One_P_y[16] = { 0x37BF51F5, 0xCBB64068, 0x6B315ECE, 0x2BCE3357, 0x7C0F9E16, 0x8EE7EB4A, 0xFE1A7F9B, 0x4FE342E2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t One_P_z[16] = { 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t Three_P_x[16] = { 0x74677826, 0xB349E309, 0x5B183A80, 0x7A08C525, 0x51FF2D3F, 0xC1639A86, 0xEF13EEB2, 0xB07647DE, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Three_P_y[16] = { 0x43238E4C, 0x81CB8368, 0x954E8ADF, 0xF69EEFDC, 0xD508DCF0, 0x1CA89A76, 0xF014119A, 0xC88A4E6D, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Three_P_z[16] = { 0x6ECDD6E2, 0xB16A0FB6, 0x4A06E794, 0x4985EC61, 0xA110D9D1, 0x9195511D, 0xABD70D36, 0x11DAA925, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t Five_P_x[16] = { 0x68DE9262, 0x0C24B048, 0x135A5F35, 0xD61EF54B, 0xE6A04259, 0x8AFCD901, 0x321B4BAC, 0x3587FE2B, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Five_P_y[16] = { 0xEFA48AA9, 0x91CE5A8C, 0x33D16005, 0x1A61E155, 0x1A0629CF, 0xB915B207, 0xAB5DFC34, 0x2BB5B613, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Five_P_z[16] = { 0x7C59A3B6, 0xF8D37FA6, 0xEE97C0ED, 0x3D8780EE, 0x4E4AB615, 0x175C41F8, 0x6577019C, 0xE242554E, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t Seven_P_x[16] = { 0x9424F2C5, 0x7F240601, 0x72BC667A, 0x2E57F932, 0xF2200721, 0x0102572E, 0xB47FE145, 0xD6FD6807, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Seven_P_y[16] = { 0xCC67CD7C, 0x664C10FC, 0x29671C3B, 0xAFC15A86, 0x8639FE12, 0x83CCBD61, 0xD27FEBC1, 0x4741C3D5, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t Seven_P_z[16] = { 0x3C7D4EF8, 0xF0F3463A, 0x6FE0509F, 0x65A0C61E, 0x2A53A02A, 0xB5A60770, 0xB7D972D6, 0x4A4F0483, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t mOne_P_x[16] = { 0xD898C296, 0xF4A13945, 0x2DEB33A0, 0x77037D81, 0x63A440F2, 0xF8BCE6E5, 0xE12C4247, 0x6B17D1F2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mOne_P_y[16] = { 0xC840AE0A, 0x3449BF97, 0x94CEA131, 0xD431CCA9, 0x83F061E9, 0x711814B5, 0x01E58065, 0xB01CBD1C, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mOne_P_z[16] = { 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t mThree_P_x[16] = { 0x74677826, 0xB349E309, 0x5B183A80, 0x7A08C525, 0x51FF2D3F, 0xC1639A86, 0xEF13EEB2, 0xB07647DE, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mThree_P_y[16] = { 0xBCDC71B3, 0x7E347C97, 0x6AB17520, 0x09611024, 0x2AF7230F, 0xE3576589, 0x0FEBEE66, 0x3775B191, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mThree_P_z[16] = { 0x6ECDD6E2, 0xB16A0FB6, 0x4A06E794, 0x4985EC61, 0xA110D9D1, 0x9195511D, 0xABD70D36, 0x11DAA925, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	uint32_t mFive_P_x[16] = { 0x68DE9262, 0x0C24B048, 0x135A5F35, 0xD61EF54B, 0xE6A04259, 0x8AFCD901, 0x321B4BAC, 0x3587FE2B, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mFive_P_y[16] = { 0x105B7556, 0x6E31A573, 0xCC2E9FFA, 0xE59E1EAB, 0xE5F9D630, 0x46EA4DF8, 0x54A203CC, 0xD44A49EB, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mFive_P_z[16] = { 0x7C59A3B6, 0xF8D37FA6, 0xEE97C0ED, 0x3D8780EE, 0x4E4AB615, 0x175C41F8, 0x6577019C, 0xE242554E, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t mSeven_P_x[16] = { 0x9424F2C5, 0x7F240601, 0x72BC667A, 0x2E57F932, 0xF2200721, 0x0102572E, 0xB47FE145, 0xD6FD6807, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mSeven_P_y[16] = { 0x33983283, 0x99B3EF03, 0xD698E3C4, 0x503EA57A, 0x79C601ED, 0x7C33429E, 0x2D80143F, 0xB8BE3C29, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t mSeven_P_z[16] = { 0x3C7D4EF8, 0xF0F3463A, 0x6FE0509F, 0x65A0C61E, 0x2A53A02A, 0xB5A60770, 0xB7D972D6, 0x4A4F0483, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };


	//4.2.1Line
	if (k_i[i] >0)//1,3,5,7
	{
		if (k_i[i] == 1)
		{
			array16_copy(One_P_x, P->X);
			array16_copy(One_P_y, P->Y);
			array16_copy(One_P_z, P->Z);
		}
		else if (k_i[i] == 3)
		{
			array16_copy(Three_P_x, P->X);
			array16_copy(Three_P_y, P->Y);
			array16_copy(Three_P_z, P->Z);
		}
		else if (k_i[i] == 5)
		{
			array16_copy(Five_P_x, P->X);
			array16_copy(Five_P_y, P->Y);
			array16_copy(Five_P_z, P->Z);
		}
		else if (k_i[i] == 7)
		{
			array16_copy(Seven_P_x, P->X);
			array16_copy(Seven_P_y, P->Y);
			array16_copy(Seven_P_z, P->Z);
		}
	}
	//4.2.2Line
	else//9, 11, 13, 15
	{
		if (k_i[i] == -7)
		{
			array16_copy(mSeven_P_x, P->X);
			array16_copy(mSeven_P_y, P->Y);
			array16_copy(mSeven_P_z, P->Z);
		}
		else if (k_i[i] == -5)
		{
			array16_copy(mFive_P_x, P->X);
			array16_copy(mFive_P_y, P->Y);
			array16_copy(mFive_P_z, P->Z);
		}
		else if (k_i[i] == -3)
		{
			array16_copy(mThree_P_x, P->X);
			array16_copy(mThree_P_y, P->Y);
			array16_copy(mThree_P_z, P->Z);
		}
		else if (k_i[i] == -1)
		{
			array16_copy(mOne_P_x, P->X);
			array16_copy(mOne_P_y, P->Y);
			array16_copy(mOne_P_z, P->Z);
		}
	}

}

//1보다 크거나 같으면 1을 return , else 0을 return
int greater_or_equal_one(uint32_t* k)
{
	int i = 0;
	for (i = 7; i > 0; i--)
	{
		if (k[i] > 0)
		{
			return 1;
		}
	}
	if (k[0] >= 1)
	{
		return 1;
	}
	return 0;
}

void Fixed_Comb_Bit_table(char (*table)[32], uint32_t *k)
{
	int i = 0, j = 0;
	
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 32; j++)
		{
			table[i][j] = (k[i] >> j) & 0x01;
		}
	}
}

void Fixed_Comb_P_table(char(*table)[32], int j, struct Affin* P)
{
	uint32_t ZERO_Line_P_X[16] = { 0xD898C296, 0xF4A13945, 0x2DEB33A0, 0x77037D81, 0x63A440F2, 0xF8BCE6E5, 0xE12C4247, 0x6B17D1F2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t ZERO_Line_P_Y[16] = { 0x37BF51F5, 0xCBB64068, 0x6B315ECE, 0x2BCE3357, 0x7C0F9E16, 0x8EE7EB4A, 0xFE1A7F9B, 0x4FE342E2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t ONE_Line_P_X[16] = { 0x185A5943, 0x3A5A9E22, 0x5C65DFB6, 0x1AB91936, 0x262C71DA, 0x21656B32, 0xAF22AF89, 0x7FE36B40, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t ONE_Line_P_Y[16] = { 0x699CA101, 0xD50D152C, 0x7B8AF212, 0x74B3D586, 0x07DCA6F1, 0x9F09F404, 0x25B63624, 0xE697D458, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t TWO_Line_P_X[16] = { 0x8E14DB63, 0x90E75CB4, 0xAD651F7E, 0x29493BAA, 0x326E25DE, 0x8492592E, 0x2811AAA5, 0x0FA822BC, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t TWO_Line_P_Y[16] = { 0x5F462EE7, 0xE4112454, 0x50FE82F5, 0x34B1A650, 0xB3DF188B, 0x6F4AD4BC, 0xF5DBA80D, 0xBFF44AE8, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t THREE_Line_P_X[16] = { 0x7512218E, 0xA84AA939, 0x74CA0141, 0xE9A521B0, 0x18A2E902, 0x57880B3A, 0x12A677A6, 0x4A5B5066, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t THREE_Line_P_Y[16] = { 0x4C4F3840, 0x0BEADA7A, 0x19E26D9D, 0x626DB154, 0xE1627D40, 0xC42604FB, 0xEAC089F1, 0xEB13461C, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t FOUR_Line_P_X[16] = { 0xD789BD85, 0x57C84FC9, 0xC297EAC3, 0xFC35FF7D, 0x88C6766E, 0xFB982FD5, 0xEEDB5E67, 0x447D739B, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t FOUR_Line_P_Y[16] = { 0x72E25B32, 0x0C7E33C9, 0xA7FAE500, 0x3D349B95, 0x3A4AAFF7, 0xE12E9D95, 0x834131EE, 0x2D4825AB, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t FIVE_Line_P_X[16] = { 0xF7F82F2A, 0xAEE9C75D, 0x4AFDF43A, 0x9E4C3587, 0x37371326, 0xF5622DF4, 0x6EC73617, 0x8A535F56, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t FIVE_Line_P_Y[16] = { 0x223094B7, 0xC5F9A0AC, 0x4C8C7669, 0xCDE53386, 0x085A92BF, 0x37E02819, 0x68B08BD7, 0x0455C084, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t SIX_Line_P_X[16] = { 0x313728BE, 0x6CF20FFB, 0xA3C6B94A, 0x96439591, 0x44315FC5, 0x2736FF83, 0xA7849276, 0xA6D39677, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t SIX_Line_P_Y[16] = { 0xC357F5F4, 0xF2BAB833, 0x2284059B, 0x824A920C, 0x2D27ECDF, 0x66B8BABD, 0x9B0B8816, 0x674F8474, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	
	uint32_t SEVEN_Line_P_X[16] = { 0xE895DF07, 0x6A703F10, 0x01876BD8, 0xFD75F3FA, 0x0CE08FFE, 0xEB5B06E7, 0x2783DFEE, 0x68F6B854, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
	uint32_t SEVEN_Line_P_Y[16] = { 0x78712655, 0x90C76F8A, 0xF310BF7F, 0xCF5293D2, 0xFDA45028, 0xFBC8044D, 0x92E40CE6, 0xCBE1FEBA, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	if (j == 0)
	{
		array16_copy(ZERO_Line_P_X, P->X);
		array16_copy(ZERO_Line_P_Y, P->Y);
	}
	else if (j == 1)
	{
		array16_copy(ONE_Line_P_X, P->X);
		array16_copy(ONE_Line_P_Y, P->Y);
	}
	else if (j == 2)
	{
		array16_copy(TWO_Line_P_X, P->X);
		array16_copy(TWO_Line_P_Y, P->Y);
	}
	else if (j == 3)
	{
		array16_copy(THREE_Line_P_X, P->X);
		array16_copy(THREE_Line_P_Y, P->Y);
	}
	else if (j == 4)
	{
		array16_copy(FOUR_Line_P_X, P->X);
		array16_copy(FOUR_Line_P_Y, P->Y);
	}
	else if (j == 5)
	{
		array16_copy(FIVE_Line_P_X, P->X);
		array16_copy(FIVE_Line_P_Y, P->Y);
	}
	else if (j == 6)
	{
		array16_copy(SIX_Line_P_X, P->X);
		array16_copy(SIX_Line_P_Y, P->Y);
	}
	else
	{
		array16_copy(SEVEN_Line_P_X, P->X);
		array16_copy(SEVEN_Line_P_Y, P->Y);
	}
}

void ECLtoR_Fixed_Comb(char(*table)[32], struct Jacobian* Q)
{
	uint32_t Z[16]= { 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };

	struct Jacobian P_Jacobian;
	struct Jacobian* Pointer_P_Jacobian = &P_Jacobian;

	struct Affin P_Affin;
	struct Affin* Pointer_P_Affin = &P_Affin;

	struct Jacobian temp_Q;
	struct Jacobian* Pointer_tQ = &temp_Q;

	
	int i = 0, j = 0, count = 0;
	
	for (i = 31; i >= 0; i--)
	{
		for (j = 0; j < 8; j++)
		{
			
			//1일때만 동작
			if (table[j][i] == 1)
			{
				//처음에는 ECADD_Jacobian 안하고 대입만 해주면된다.
				if (count == 0)
				{
					Fixed_Comb_P_table(table, j, Pointer_P_Affin);
					
					//Jacobian Q로
					array16_copy(Pointer_P_Affin->X, Q->X);
					array16_copy(Pointer_P_Affin->Y, Q->Y);
					array16_copy(Z, Q->Z);

					count = 1;
				}
				//두번째부터는 ECadd_Jacobian 해주어야한다.
				else
				{
					Fixed_Comb_P_table(table, j, Pointer_P_Affin);

					array16_copy(Q->X, Pointer_tQ->X);
					array16_copy(Q->Y, Pointer_tQ->Y);
					array16_copy(Q->Z, Pointer_tQ->Z);

					ECADD_Jacobian(Pointer_tQ, Pointer_P_Affin, Q);
				}
			}
		}
		//하나의 열 ECADD끝난후 ECDBL_Jacobian 해주면된다.
		if ((count == 1)&&(i!=0))
		{
			array16_copy(Q->X, Pointer_tQ->X);
			array16_copy(Q->Y, Pointer_tQ->Y);
			array16_copy(Q->Z, Pointer_tQ->Z);
			ECDBL_Jacobian(Pointer_tQ, Q);
		}
	}

	//Q를 Affin으로 바꾸어주어야한다.
	change_affin(Q);
}
