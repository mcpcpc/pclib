#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	int m;
	int n;
	double ** v;
} mat_t, *mat;

mat matrixAllocate(int m, int n) {
	mat x = malloc(sizeof(mat_t));
	x->v = malloc(sizeof(double*) * m);
	x->v[0] = calloc(sizeof(double), m * n);
	for (int i = 0; i < m; i++) {
		x->v[i] = x->v[0] + n * i;
	}
	x->m = m;
	x->n = n;
	return x;
}

void matrixFree(mat m) {
	free(m->v[0]);
	free(m->v);
	free(m);
}

double matrixMean(mat m, int dim) {
	double sum = 0;
	for (int i = 0; i < m->m; i++) {
		sum += m->v[i][dim];
	}
	double x = sum / m->m;
	return x;
}

double matrixStdDevP(mat m, int dim) {
	double sum = 0;
	double mean = matrixMean(m, dim);
	for (int i = 0; i < m->m; i++) {
		sum += pow(m->v[i][dim] - mean, 2);
	}
	double sd = sqrt(sum / (m->m - 1));
	return sd;
}

mat matrixStandardize(mat x) {
	mat m = matrixAllocate(x->m, x->n);
	double mean = 0;
	double sdev = 0;
	for (int i = 0; i < m->m; i++) {
		for (int j = 0; j < m->n; j++) {
			mean = matrixMean(x, j);
			sdev = matrixStdDevP(x, j);
			m->v[i][j] = (x->v[i][j] - mean) / sdev;
		}
	}
	return m;
}

double matrixCovariance(mat m, int x, int y) {
	double sum = 0;
	double meanX = matrixMean(m, x);
	double meanY = matrixMean(m, y);
	for (int i = 0; i < m->m; i++) {
		sum += (m->v[i][x] - meanX) * (m->v[i][y] - meanY);
	}
	double cov = sum / m->m;
	return cov;
}

mat covarianceMatrix(mat x) {
	mat m = matrixAllocate(x->n, x->n);
	for (int i = 0; i < m->n; i++) {
		for (int j = i; j < m->n; j++) {
			m->v[i][j] = matrixCovariance(x, i, j);
			m->v[j][i] = m->v[i][j];
		}
	}
	return m;
}

void matrixTranspose(mat m) {
	for (int i = 0; i < m->m; i++) {
		for (int j = 0; j < i; j++) {
			double t = m->v[i][j];
			m->v[i][j] = m->v[j][i];
			m->v[j][i] = t;
		}
	}
}

mat matrixCopy(int n, double a[][n], int m) {
	mat x = matrixAllocate(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			x->v[i][j] = a[i][j];
		}
	}
	return x;
}

mat matrixMultiply(mat x, mat y) {
	if (x->n != y->m) return 0;
	mat r = matrixAllocate(x->m, y->n);
	for (int i = 0; i < x->m; i++) {
		for (int j = 0; j < y->n; j++) {
			for (int k = 0; k < x->n; k++) {
				r->v[i][j] += x->v[i][k] * y->v[k][j];
			}
		}
	}
	return r;
}

mat matrixMinor(mat x, int d) {
	mat m = matrixAllocate(x->m, x->n);
	for (int i = 0; i < d; i++) {
		m->v[i][i] = 1;
	}
	for (int i = d; i < x->m; i++) {
		for (int j = d; j < x->n; j++) {
			m->v[i][j] = x->v[i][j];
		}
	}
	return m;
}

/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], int n) {
	for (int i = 0; i < n; i++) {
		c[i] = a[i] + s * b[i];
	}
	return c;
}

/* m = I - v v^T */
mat vectorMultiply(double v[], int n) {
	mat x = matrixAllocate(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x->v[i][j] = -2 *  v[i] * v[j];
		}
	}
	for (int i = 0; i < n; i++) {
		x->v[i][i] += 1;
	}
	return x;
}

/* ||x|| */
double vnorm(double x[], int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += x[i] * x[i];
	}
	return sqrt(sum);
}

/* y = x / d */
double *vectorDivide(double x[], double d, double y[], int n) {
	for (int i = 0; i < n; i++) {
		y[i] = x[i] / d;
	}
	return y;
}

/* take c-th column of m, put in v */
double *mcol(mat m, double *v, int c) {
	for (int i = 0; i < m->m; i++) {
		v[i] = m->v[i][c];
	}
	return v;
}

void matrixPrint(mat m) {
	for(int i = 0; i < m->m; i++) {
		for (int j = 0; j < m->n; j++) {
			printf(" %8.3f", m->v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matrixHouseholder(mat m, mat *R, mat *Q) {
	mat q[m->m];
	mat z = m, z1;
	for (int k = 0; k < m->n && k < m->m - 1; k++) {
		double e[m->m], x[m->m], a;
		z1 = matrixMinor(z, k);
		if (z != m) {
			matrixFree(z);
		}
		z = z1;
		mcol(z, x, k);
		a = vnorm(x, m->m);
		if (m->v[k][k] > 0) {
			a = -a;
		}
		for (int i = 0; i < m->m; i++) {
			e[i] = (i == k) ? 1 : 0;
		}
		vmadd(x, e, a, e, m->m);
		vectorDivide(e, vnorm(e, m->m), e, m->m);
		q[k] = vectorMultiply(e, m->m);
		z1 = matrixMultiply(q[k], z);
		if (z != m) {
			matrixFree(z);
		}
		z = z1;
	}
	matrixFree(z);
	*Q = q[0];
	*R = matrixMultiply(q[0], m);
	for (int i = 1; i < m->n && i < m->m - 1; i++) {
		z1 = matrixMultiply(q[i], *Q);
		if (i > 1) {
			matrixFree(*Q);
		}
		*Q = z1;
		matrixFree(q[i]);
	}
	matrixFree(q[0]);
	z = matrixMultiply(*Q, m);
	matrixFree(*R);
	*R = z;
	matrixTranspose(*Q);
}

/* matrix of diagonal ones */
mat matrixEye(int m) {
	mat x = matrixAllocate(m, m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (i == j) {
				x->v[i][j] = 1;
			} else {
				x->v[i][j] = 0;
			}
		}
	}
	return x;
}

/* Ak+1 = Rk * Qk */
mat eigenValues(mat m, int k) {
	mat Q, R;
	matrixHouseholder(m, &R, &Q);
	mat A = matrixMultiply(R, Q);
	mat pQ = matrixEye(m->m);
	pQ = matrixMultiply(pQ, Q);
	for (int i = 1; i < k; i++) {
		matrixHouseholder(A, &R, &Q);
		pQ = matrixMultiply(pQ, Q);
		A = matrixMultiply(R, Q);
	}
	matrixFree(Q);
	matrixFree(R);
	matrixFree(A);
	return pQ;
}

mat matrixExtract(mat m, int x1, int x2, int y1, int y2) {
	int mNew = x2 - x1;
	int nNew = y2 - y1;
	mat x = matrixAllocate(mNew, nNew);
	for (int i = 0; i < mNew; i++) {
		for (int j = 0; j < nNew; j++) {
			x->v[i][j] = m->v[x1 + i][y1 + j];
		}
	}
	return x;
}

int n_components = 2;
double src[][4] = {{1,2,3,4},{5,5,6,7},{1,4,2,3},{5,3,2,1},{8,1,2,2}};

int main() {
	/* covariance matrix */
	mat x = matrixCopy(4,src, 5);
	puts("IN");
	matrixPrint(x);
	mat b = matrixStandardize(x);
	mat c = covarianceMatrix(b);
	mat eigen = eigenValues(c, 100);
	mat s = matrixExtract(eigen, 0, eigen->m, 0, n_components);
	mat pc = matrixMultiply(b, s);
	puts("PRINCIPAL COMPONENTS");
	matrixPrint(pc);

	matrixFree(s);
	matrixFree(pc);
	matrixFree(eigen);
	matrixFree(x);
	matrixFree(b);
	matrixFree(c);
	return 0;
}
