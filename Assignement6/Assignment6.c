#include <stdio.h>
#include <time.h>
#include <math.h>
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// uniform distribution
float uniform(long *idum)
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

// gaussian distribution
float gasdev(long *idum)
{
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;

	if (*idum < 0) iset = 0;
	if (iset == 0) {
		do {
			v1 = 2.0*uniform(idum) - 1.0;
			v2 = 2.0*uniform(idum) - 1.0;
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
	int i, GEN;
	scanf("%d", &GEN);
	long t = -(long)time(NULL);
	FILE* fd = fopen("result.txt", "w");
    // uniform distribution을 이용한 랜덤넘버
	for (i = 0; i < GEN; i++) {
		fprintf(fd, "%f\n", uniform(&t) * 5.0f - 3.0f);
	}
	fclose(fd);
	fd = fopen("gaussian.txt", "w");
    // gaussian distribution을 이용한 랜덤넘버
	for (i = 0; i < GEN; i++) {
        // 3s초과
        if (1.5f * gasdev(&t) + 0.5f > 1.5f*3){
            fprintf(fd, "%f\n", 4.5f);
        }
        // 3s이하
        else if (1.5f * gasdev(&t) + 0.5f < -1.5f*3){
            fprintf(fd, "%f\n", -4.5f);
        }
        else{
            fprintf(fd, "%f\n", 1.5f * gasdev(&t) + 0.5f);
        }
	}
	fclose(fd);
	return 0;
}
