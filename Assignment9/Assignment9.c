// written by odysseyj
// HYU_2015004239 정성운
// https://github.com/OdysseyJ

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*
#define TINY 1.0e-20

extern int err = 0;
extern float det = 1.0;

void nrerror(char msg[]);
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void ludcmp(double a[][4], int n, int *indx, float *d);
void lubksb(double a[][4], int n, int *indx, double b[]);

int main(void) {
	double x[200], y[200], x1[200], y1[200];

    double J_tranJ[4][4], J[200][4], tranJ[4][200];
    
    double sol[4];
    char fileName[50];
    int indx[4];
    int count = 0;
    int i, j, k;
	double temp;
	float d;
    
    // 파일 읽기
    printf("filename: ");
    scanf("%s", fileName);
    FILE *fp = fopen(fileName, "r");
    while (fscanf(fp, "%lf", &temp) != EOF) x[count] = temp, fscanf(fp, "%lf %lf %lf", &y[count], &x1[count], &y1[count]), count++;
    // 파일 읽기 종료
    
    // Jacovian Transpose, Jacovian matrix 채우기
    for (i = 0; i < count; i++) {
        tranJ[1][i + 1] = J[i + 1][1] = x[i];
        tranJ[2][i + 1] = J[i + 1][2] = y[i];
        tranJ[3][i + 1] = J[i + 1][3] = 1;
    }
    
    // J^T * J  채우기
    for (i = 1; i <= 3; i++){
        for (j = 1; j <= 3; j++){
            J_tranJ[i][j] = 0;
        }
    }
    for (i = 1; i <= 3; i++){
        for (j = 1; j <= 3; j++){
            for (k = 1; k <= count; k++){
                J_tranJ[i][j] += tranJ[i][k] * J[k][j];
            }
        }
    }
    
    // 우항 만들기.
	sol[1] = sol[2] = sol[3] = 0;
	for (i = 1; i <= count; i++) {
		sol[1] += x[i - 1] * x1[i - 1];
		sol[2] += y[i - 1] * x1[i - 1];
		sol[3] += x1[i - 1];
	}
    
    // J_tranJ 매트릭스 다시 분해하기. (LU backsubstitution 사용을 위해)
	ludcmp(J_tranJ, 3, indx, &d);
    
    // lu back substitution으로 정답 구하기
	lubksb(J_tranJ, 3, indx, sol);
    
    printf("Linear Data Fitting\n");
    printf("Solutions)\n");
    for (i = 1; i <= 3; i++) {
        printf("a[%d]: %lf\n", i, sol[i]);
    }
    
    // 우항 만들기.
	sol[1] = sol[2] = sol[3] = 0;
	for (j = 1; j <= count; j++) {
		sol[1] += x[j - 1] * y1[j - 1];
		sol[2] += y[j - 1] * y1[j - 1];
		sol[3] += y1[j - 1];
	}
    
    // J_tranJ 매트릭스 다시 분해하기.
	lubksb(J_tranJ, 3, indx, sol);
    
    // LU back substitution으로 정답 구하기.
    for (i = 1; i <= 3; i++) {
        printf("a[%d]: %lf\n", i + 3, sol[i]);
    }
	return 0;
}

void nrerror(char msg[]) {
    printf("%s", msg);
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

// LU_decomposition
void ludcmp(double a[][4], int n, int *indx, float *d)
{
    int i, imax, j, k;
    float big, dum, sum, temp;
    float *vv;
    vv = vector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big) big = temp;
        if (big == 0.0) {
            printf("Singular matrix in routine ludcmp\n");
            err++;
            return;
        }
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i<j; i++) {
            sum = a[i][j];
            for (k = 1; k<i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k<j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++) a[i][j] *= dum;
        }
    }
    free_vector(vv, 1, n);
}

// LU_backsubstitution
void lubksb(double a[][4], int n, int *indx, double b[])
{
    int i, ii = 0, ip, j;
    float sum;
    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}
