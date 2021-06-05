#include <unistd.h>
#include <math.h>
#include <stdio.h>

double mean(size_t m, size_t n, double mat[][n], size_t d) {
	double sum = 0;
	for (size_t i = 0; i < m; i++) {
		sum += mat[i][d];
	}
	double xBar = sum / m;
	return xBar;
}

double stdDev(size_t m, size_t n, double mat[][n], size_t d) {
	double sum = 0;
	double xBar = mean(m, n, mat, d);
	for (size_t i = 0; i < m; i++) {
		sum += pow(mat[i][d] - xBar, 2);
	}
	double sd = sqrt(sum / (m - 1));
	return sd;
}


double covariance(size_t m, size_t n, double mat[][n], size_t x, size_t y) {
	double sum = 0;
	double meanX = mean(m, n, mat, x);
	double meanY = mean(m, n, mat, y);
	for (size_t i = 0; i < m; i++) {
		sum += (mat[i][x] - meanX) * (mat[i][y] - meanY);
	}
	double cov = sum / m;
	return cov;
}

void covarianceMatrix(size_t m, size_t n, double src[][n], double dest[][n]) {
	for (size_t x = 0; x < n; x++) {
		for (size_t y = x; y < n; y++) {
			dest[x][y] = covariance(m, n, src, x, y);
			dest[y][x] = dest[x][y];
		}
	}
}

void printMat(size_t m, size_t n, double mat[][n]) {
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			printf("\t%f", mat[i][j]);
		}
		puts("\n");
	}
}

void standardMatrix(size_t m, size_t n, double src[][n], double dest[][n]) {
	double meanCol = 0;
	double stdDCol = 0;
	for (size_t x = 0; x < m; x++) {
		for (size_t y = 0; y < n; y++) {
			meanCol = mean(m, n, src, y);
			stdDCol = stdDev(m, n, src, y);
			dest[x][y] = (src[x][y] - meanCol) / stdDCol;
		}
	}
}
/*
void copyMatrixCol(size_t, size_m, size_t n, double src[][n], double dest[], size_t i) {

}

void qrDecompose(size_t m, size_t n, double src[][n], double dest[][n]) {
	 double T[m];
	 double S[m];
	 for (size_t y = 0; y < n; y++) {
		for (size_t x = 0; x < y; x++) {
		}
	 }
}
*/
int main() {
	// assumption is that cols (n) is less than rows (m)
	size_t dims[2] = {5, 4};
	double src[5][4] = {{1,2,3,4},{5,5,6,7},{1,4,2,3},{5,3,2,1},{8,1,2,2}};
	double std[dims[0]][dims[1]];
	standardMatrix(dims[0], dims[1], src, std);
	size_t newDim = 4;
	double dest[newDim][newDim];
	covarianceMatrix(dims[0], dims[1], std, dest);
	puts("covariance matrix:");
	printMat(newDim, newDim, dest);
	puts("QR algorithm:");

}
