//const char* dgemm_desc = "My awesome matmul.";

#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#ifndef L1_BS
#define L1_BS ((int) 16*16)
#endif

#ifndef L2_BS
#define L2_BS ((int) 8)
#endif

// This code was written by Marc Aurele Gilles for the Matrix multiply hw for cs5220 at Cornell



/*

  A,B,C are M-by-M
  L1_BS is the size of sub matrix that will fit into L1
  L2_BS is the number of submatrix that will fit in L2
  L3*L2 is the number of submatrix that will fit in L3

*/
void row_to_block(const int M, const int padsize,  const int nblock, const int *restrict A, int *restrict newA)
{
	// converts to block indexing and pads the new matrix with zeros so that it is divisble by L1_BS
	int bi,bj,i,j;
	int inf = padsize+1;
	for(bi=0; bi < nblock; ++bi){
		for(bj=0; bj < nblock; ++bj){
			for(i=0; i < L1_BS; ++i){
				for(j=0; j < L1_BS; ++j){
					if ((bj*L1_BS+j) >= M || bi*L1_BS+i >= M){ 
						// we can optimize this to delete this "if". this is doing the padding 
						newA[((bj*nblock+bi)*L1_BS*L1_BS+ i*L1_BS+j)]=inf;
					}
					else{
						newA[(bj*nblock+bi)*L1_BS*L1_BS+ i*L1_BS+j]=A[(j+bj*L1_BS)*M+bi*L1_BS+i]; 
					}				
				}
			}
		}
	}
}








void block_to_row(const int M, const int nblock, int *restrict A, const int *restrict newA)
{
	int bi, bj,i,j;	
	for(bi=0; bi < nblock; ++bi){
		for(bj=0; bj < nblock; ++bj){
			for(i=0; i < L1_BS; ++i){
				for(j=0; j < L1_BS; ++j){
					if ((bj*L1_BS+j)>= M || bi*L1_BS+i >= M){
					   	// we can optimize this to delete this "if". this is doing the padding 
					}
					else{
						A[(j+bj*L1_BS)*M+bi*L1_BS+i]=newA[(bj*nblock+bi)*L1_BS*L1_BS+ i*L1_BS+j];
				   	}

				}
			}
		}
	}
}
/* not used in current
void row_to_block_transpose(const int M, const int nblock, const int *restrict A, int *restrict newA)
{
	// converts to block indexing and pads the new matrix with zeros so that it is divisble by L1_BS
	int bi,bj,i,j;	
	for(bi=0; bi < nblock; ++bi){
		for(bj=0; bj < nblock; ++bj){
			for(i=0; i < L1_BS; ++i){
				for(j=0; j < L1_BS; ++j){
					if ((bj*L1_BS+j) >= M || bi*L1_BS+i >= M){ 
						// we can optimize this to delete this "if". this is doing the padding 
						newA[((bj*nblock+bi)*L1_BS*L1_BS+ j*L1_BS+i)]=0;
					}
					else{
						newA[(bj*nblock+bi)*L1_BS*L1_BS+ j*L1_BS+i]=A[(j+bj*L1_BS)*M+bi*L1_BS+i]; 
					}				
				}
			}
		}
	}
}
*/






	
int do_block(const int M, const int nblock,
              const int * restrict A, int * restrict C,
              const int bi, const int bj, const int bk)
{	// A is old matrix, C is new matrix
	int i, j, k, BA_Ar, BA_Ac, sub_BA_Ar, sub_BA_Ac, BA_C, sub_BA_C;
    __assume_aligned(A, 64);
    __assume_aligned(C, 64);
   int done =1;
   printf("nblock: %d", nblock);
	// BA stands for block adress 
    BA_Ar=(bk*nblock+bi)*L1_BS*L1_BS; // for the column block of old matrix
    BA_Ac=(bj*nblock+bk)*L1_BS*L1_BS; // for the row block of old matrix
    BA_C=(bj*nblock+bi)*L1_BS*L1_BS;  // for the new matrix
    for (i = 0; i < L1_BS; ++i) {
	// finds sub_BA, tells compiler its aligned
	sub_BA_Ar=BA_Ar+L1_BS*i;
	__assume(sub_BA_Ar%128==0);
	//same for C                
 	sub_BA_C=BA_C+L1_BS*i;
	__assume(sub_BA_C%128==0);

	for (j = 0; j < L1_BS; ++j){
	    int cij = C[sub_BA_C+j];
            for (k = 0; k < L1_BS; ++k) {
                __assume(sub_BA_Ac%128==0); // need to change factor since int are smaller than float
                sub_BA_Ac=BA_Ac+L1_BS*k;
		if ( sub_BA_Ar+ k >= nblock*L1_BS*nblock*L1_BS ){ printf("BA_Ar %d, BA_Ac %d, BA_C %d, bi : %d, bj: %d, bk %d,  sub_BA_Ar+ k: %d, pad size:%d \n "
,BA_Ar, BA_Ac, BA_C, bi, bj, bk, sub_BA_Ar +k, nblock*L1_BS*nblock*L1_BS);
		assert (sub_BA_Ar+ k< nblock*L1_BS);  }
		assert( sub_BA_Ac+j < nblock*L1_BS*nblock*L1_BS); 
                if (A[sub_BA_Ar+k] + A[sub_BA_Ac+j] < cij){
			cij = A[sub_BA_Ar+k] + A[sub_BA_Ac+j];
                   	 done=0;
				}
		}
	 C[sub_BA_C+j]= cij;
        }
    }
}


void setup_indices(
	const int M, 
	const int* A, const int* C, 
	int* bA, int* bC, 
	int* pad_size, int* nblock, 
	int* L1nblock, int* L2nblock, 
	int* rem
) {
	/* pad size : size of matrix after padding
	   block i= row index of block
	   block j = column index of block
	   nblock = number of blocks
	   L2bi, L2bj index of L2
           L2nblock= number of L1block that fit into L2
	   L2rem number of remaining blocks 
           rem : remainder of blocks after l2 blocking
	   Ai, Aj, Ak : "addresses of i j k loops"
	*/
	printf("starting setup \n");
	printf(" M: %d, L1_BS:%d M/L1_BS: %d\n ", M, L1_BS, M/L1_BS);
	if (M%L1_BS==0){
		*nblock=M/L1_BS;
		*pad_size=M;
	}
	else{ 
		*pad_size=((M/L1_BS)+1)*L1_BS;
		*nblock=M/L1_BS+1;
	}
	
	// number of L2
	if(*pad_size%(L2_BS*L1_BS)==0){
	*L2nblock=*pad_size/(L2_BS*L1_BS);
	// define remainder to be a whole block in this case for consistency
	 *rem = L2_BS;}

	else{*L2nblock=*pad_size/(L2_BS*L1_BS)+1;
	     *rem= (*pad_size%(L2_BS*L1_BS))/L1_BS;	}


	printf("done allocating \n");
	bA= (int*) _mm_malloc((*pad_size)*(*pad_size)*sizeof(int),64);
	bC= (int*) _mm_malloc((*pad_size)*(*pad_size)*sizeof(int),64);
	printf("padsize: %d \n", *pad_size);
	// change indexing
	row_to_block(M, *nblock, *pad_size,  A, bA);
	row_to_block(M, *nblock, *pad_size, C, bC);
    	
  }

/* 
 *
 */
int square_dgemm(const int M,
		 const int *bA, int *bC,
		 const int  pad_size,const int nblock,const  int L1nblock, const int L2nblock, const int rem)
    {

    int bi, bj, bk, L2bi, L2bj, L2bk;
    int done=1;
	// MAIN LOOP
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bk=0; L2bk < L2nblock-1; ++L2bk){
	for (L2bj=0; L2bj < L2nblock-1; ++L2bj){
	for (L2bi=0; L2bi < L2nblock-1; ++L2bi){
	for (bk = 0; bk < L2_BS; ++bk) {
		int Ak=L2bk*L2_BS +bk;
		for (bj = 0; bj < L2_BS; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < L2_BS; ++bi) {
				int Ai=L2bi*L2_BS+bi;
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
			}
			}
		}
	}
	}
	}

	printf("finished main loop \n");

		

	// ADDITIONAL LOOPS TO AVOID IF STATEMENTS
	// there are 8 cases: (we denote j not at boundary by j0 and j at boundary with j1
	// 1: k0 j0 i0 (main loop)
	// 2: k1 j0 i0
	// 3: k0 j1 i0
	// 4: k0 j0 i1
	// 5: k1 j1 i0
	// 6: k1 j0 i1
	// 7: k0 j0 i1
	// 8: k1 j1 i1
	//printf"we got here\n");
	

	// case 2 k1 j0 i0
	L2bk=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bj=0; L2bj < L2nblock-1; ++L2bj){
	for (L2bi=0; L2bi < L2nblock-1; ++L2bi){
	for (bk = 0; bk < rem; ++bk) {
		int Ak=L2bk*L2_BS+bk;
		for (bj = 0; bj < L2_BS; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < L2_BS; ++bi) {
				int Ai=L2bi*L2_BS+bi;	
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
           
			
		}
	}
	}
	}
	}
    	printf("done with case 2\n");

	// case 3: k0 j1 i0
	L2bj=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bk=0; L2bk < L2nblock-1; ++L2bk){
	for (L2bi=0; L2bi < L2nblock-1; ++L2bi){
	for (bk = 0; bk < L2_BS; ++bk) {
		int Ak=L2bk*L2_BS+bk;
		for (bj = 0; bj < rem; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < L2_BS; ++bi) {
				int Ai=L2bi*L2_BS+bi;	
				done=done&&do_block(M, nblock, bA,  bC, Ai, Aj, Ak);
                                
			
		}
	}
	}
	}
	}
	printf("done with case 3");
	// case 4: k0 j0 i1
	L2bi=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bk=0; L2bk < L2nblock-1; ++L2bk){
	for (L2bj=0; L2bj < L2nblock-1; ++L2bj){
	for (bk = 0; bk < L2_BS; ++bk) {
		int Ak=L2bk*L2_BS+bk;
		for (bj = 0; bj < L2_BS; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < rem; ++bi) {
				int Ai=L2bi*L2_BS+bi;	
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
                               
			
		}
	}
	}
	}
	}
	printf("done with case 4");
	

	// case 5: k1 j1 i0
	L2bj=L2nblock-1;
	L2bk=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bi=0; L2bi < L2nblock-1; ++L2bi){
	for (bk = 0; bk < rem; ++bk) {
		int Ak=L2bk*L2_BS +bk;
		for (bj = 0; bj < rem; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < L2_BS; ++bi) {

				int Ai=L2bi*L2_BS+bi;	               
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
			}		
		}
	}
	}
	
	// case 6: k1 j0 i1
	L2bi=L2nblock-1;
	L2bk=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bj=0; L2bj < L2nblock-1; ++L2bj){
	for (bk = 0; bk < rem; ++bk) {
		int Ak=L2bk*L2_BS +bk;
		for (bj = 0; bj < L2_BS; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < rem; ++bi) {

				int Ai=L2bi*L2_BS+bi;	
                       
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
			}		
		}
	}
	}
	// case 7: k0 j1 i1
	L2bj=L2nblock-1;
	L2bi=L2nblock-1;
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (L2bk=0; L2bk < L2nblock-1; ++L2bk){
	for (bk = 0; bk < L2_BS; ++bk) {
		int Ak=L2bk*L2_BS +bk;
		for (bj = 0; bj < rem; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < rem; ++bi) {

				int Ai=L2bi*L2_BS+bi;	
                       
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
			}		
		}
	}
	}
	
	// case 8: k1 j1 i1
	L2bk=L2nblock-1;
	L2bj=L2nblock-1;
 	L2bi=L2nblock-1;
	printf("Ai = %d, pad_size %d \n", L2bk*L2_BS + rem-1, pad_size);
    #pragma omp parallel for shared(bA, bC) reduction(&& : done)
	for (bk = 0; bk < rem; ++bk) {
		int Ak=L2bk*L2_BS +bk;
		for (bj = 0; bj < rem; ++bj) {
			int Aj=L2bj*L2_BS+bj;
			for (bi = 0; bi < rem; ++bi) {
				int Ai=L2bi*L2_BS+bi;	

        
				done=done&&do_block(M, nblock, bA, bC, Ai, Aj, Ak);
			
		}
	}
	}
	
	
	assert(0);	
	return done;
}

