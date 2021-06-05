#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
typedef struct {
	size_t m
	size_t n;
	double **v;
} matrix_t, *matrix;
 
matrix matrixAdd(size_t m, size_t n) {
	matrix mat = malloc(sizeof(matrix_t));
	mat->v = malloc(sizeof(double*) * m);
	mat->v[0] = calloc(sizeof(double), m * n);
	for (size_t i = 0; i < m; i++) {
		mat->v[i] = mat->v[0] + n * i;
	}
	mat->m = m;
	mat->n = n;
	return mat;
}

void matrixDelete(matrix mat) {
	free(mat->v[0]);
	free(mat->v);
	free(mat);
}

void matrixTranspose(matrix mat) {
	for (size_t i = 0; i < mat->m; i++) {
		for (size_t j = 0; j < i; j++) {
			double t = mat->v[i][j];
			mat->v[i][j] = mat->v[j][i];
			mat->v[j][i] = t;
		}
	}
}

matrix matrixCopy(size_t m, size_t n, double a[][n]) {
	matrix ret = matrixNew(m, n);
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			ret->v[i][j] = a[i][j];
		}
	}
	return ret;
}

matrix matrixMultiply(matrix mat1, matrix mat2) {
	if (mat1->n != mat2->m) {
		return 0;
	}
	matrix ret = matrixNew(mat1->m, mat2->n);
	for (size_t i = 0; i < mat1->m; i++)
		for (size_t j = 0; j < mat2->n; j++)
			for (size_t k = 0; k < mat1->n; k++)
				ret->v[i][j] += mat1->v[i][k] * mat2->v[k][j];
	return ret;
}

matrix matrixMinor(matrix mat, size_t d) {
	matrix ret = matrixNew(mat->m, mat->n);
	for (size_t i = 0; i < d; i++) {
		ret->v[i][i] = 1;
	}
	for (size_t i = d; i < mat->m; i++) {
		for (size_t j = d; j < mat->n; j++) {
			ret->v[i][j] = mat->v[i][j];
	}
	return ret;
}

/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], size_t n) {
	for (size_t i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
}

/* m = I - v v^T */
matrix vmul(double v[], size_t n)
{
	matrix x = matrix_new(n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			x->v[i][j] = -2 *  v[i] * v[j];
		}
	}
	for (size_t i = 0; i < n; i++) {
		x->v[i][i] += 1;
	}
	return x;
}

/* ||x|| */
double vnorm(double x[], size_t n) {
	double sum = 0;
	for (size_t i = 0; i < n; i++) {
		sum += x[i] * x[i];
	}
	return sqrt(sum);
}

double* vdiv(double x[], double d, double y[], size_t n) {
	for (size_t i = 0; i < n; i++) {
		y[i] = x[i] / d;
	}
	return y;
}

/* take c-th column of mat, put in v */
double* mcol(matrix mat, double *v, int c) {
	for (size_t i = 0; i < mat->m; i++)
		v[i] = mat->v[i][c];
	return v;
}

void matrixPrint(matrix m) {
	for(size_t i = 0; i < mat->m; i++) {
		for (size_t j = 0; j < mat->n; j++) {
			printf(" %8.3f", mat->v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void householderTranform(matrix mat, matrix *R, matrix *Q) {
	matrix q[mat->m];
	matrix z = mat, z1;
	for (size_t k = 0; (k < mat->n) && k < (mat->m - 1); k++) {
		double e[mat->m], x[mat->m], a;
		z1 = matrixMinor(z, k);
		if (z != mat) {
			matrixDelete(z);
		}
		z = z1;
		mcol(z, x, k);
		a = vnorm(x, mat->m);
		if (mat->v[k][k] > 0) {
			a = -a;
 		}
		for (size_t i = 0; i < mat->m; i++) {
			e[i] = (i == k) ? 1 : 0;
 		}
		vmadd(x, e, a, e, mat->m);
		vdiv(e, vnorm(e, mat->m), e, mat->m);
		q[k] = vmul(e, mat->m);
		z1 = matrixMultiply(q[k], z);
		if (z != mat) {
			matrixDelete(z);
		}
		z = z1;
	}
	matrixDelete(z);
	*Q = q[0];
	*R = matrixMultiply(q[0], mat);
	for (size_t i = 1; i < mat->n && i < mat->m - 1; i++) {
		z1 = matrixMultiply(q[i], *Q);
		if (i > 1) {
			matrixDelete(*Q);
		}
		*Q = z1;
		matrixDelete(q[i]);
	}
	matrixDelete(q[0]);
	z = matrixMultiply(*Q, mat);
	matrixDelete(*R);
	*R = z;
	matrixTranspose(*Q);
}

double in[][3] = {
	{ 12, -51,   4},
	{  6, 167, -68},
	{ -4,  24, -41},
	{ -1, 1, 0},
	{ 2, 0, 3},
};

int main() {
	matrix R, Q;
	matrix x = matrixCopy(3, in, 5);
	householder(x, &R, &Q);
 
	puts("Q"); matrixShow(Q);
	puts("R"); matrixShow(R);
 
	// to show their product is the input matrix
	matrix mat = matrixMultiply(Q, R);
	puts("Q * R"); matrixShow(mat);
 
	matrixDelete(x);
	matrixDelete(R);
	matrixDelete(Q);
	matrixDelete(mat);
	return 0;
}
