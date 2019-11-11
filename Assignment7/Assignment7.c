#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <stdlib.h>
#define SIZE 11
#define NR_END 1
#define FREE_ARG char*
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	printf("%s", error_text);
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void jacobi(float a[][SIZE + 1], int n, float d[], float v[][SIZE + 1], int *nrot)
{
	int j, iq, ip, i;
	float tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

	b = vector(1, n);
	z = vector(1, n);
	for (ip = 1; ip <= n; ip++) {
		for (iq = 1; iq <= n; iq++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip = 1; ip <= n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	*nrot = 0;
	for (i = 1; i <= 50; i++) {
		sm = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z, 1, n);
			free_vector(b, 1, n);
			return;
		}
		if (i < 4)
			tresh = 0.2*sm / (n*n);
		else
			tresh = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++) {
				g = 100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip]) + g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq]) + g) == (float)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ((float)(fabs(h) + g) == (float)fabs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t*t);
					s = t*c;
					tau = s / (1.0 + c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j = 1; j <= ip - 1; j++) {
						ROTATE(a, j, ip, j, iq)
					}
					for (j = ip + 1; j <= iq - 1; j++) {
						ROTATE(a, ip, j, j, iq)
					}
					for (j = iq + 1; j <= n; j++) {
						ROTATE(a, ip, j, iq, j)
					}
					for (j = 1; j <= n; j++) {
						ROTATE(v, j, ip, j, iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip = 1; ip <= n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE

void eigsrt(float d[], float v[][SIZE + 1], int n)
{
	int k, j, i;
	float p;

	for (i = 1; i<n; i++) {
		p = d[k = i];
		for (j = i + 1; j <= n; j++)
			if (d[j] >= p) p = d[k = j];
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (j = 1; j <= n; j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran1(long *idum)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA*(*idum - k*IQ) - IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if (*idum < 0) *idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
float gasdev(long *idum)
{
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	if (*idum < 0) iset = 0;
	if (iset == 0) {
		do {
			v1 = 2.0*ran1(idum) - 1.0;
			v2 = 2.0*ran1(idum) - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq) / rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else {
		iset = 0;
		return gset;
	}
}

int main(void) {
	int i, j;
	float m[SIZE + 1][SIZE + 1];
	float v[SIZE + 1][SIZE + 1];
	int nrot;
	float eigen[SIZE + 1];
    float arr[SIZE+1];
	long t = -(long)time(NULL);
    // 랜덤 넘버 매트릭스 만들기 (symmetric)
	for (i = 1; i <= SIZE; i++) for (j = i; j <= SIZE; j++) m[i][j] = m[j][i] = gasdev(&t);
    
    // jacobi funciton
	jacobi(m, SIZE, eigen, v, &nrot);
    
    // 소팅하기.
	eigsrt(eigen, v, SIZE);
    
    // 출력.
	for (i = 1; i <= SIZE; i++) {
        printf("-------------%d------------\n",i);
		printf("eigenvalue: \n%f\n", eigen[i]);
        printf("eigenvector in the descending order:\n");
        for (j = 1; j <= SIZE; j++) {
            printf("%f\n", v[j][i]);
        }
	}
	return 0;
}
