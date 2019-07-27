#pragma once

#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

// Gauss Laguerre abscissas x and weights w
// \int_0^\infty x^alpha e^{-x} f(x) dx = \sum_0^n w_i f(x_i)
// precision test : see  https://keisan.casio.com/exec/system/1281279441

class GHGaussL
{

public:
	GHGaussL();
	~GHGaussL();

	inline static void getList(long n,  double alpha,  std::vector<double>& x, std::vector<double>& w) {
		const long MAXIT = 10;
		const double EPS = 1.0e-31;
		const double MIN = 1.0e-199;
		long i, its, j;
		double ai, p1, p2, p3, pp, z, z1;

		for (i = 0; i < n; i++) {
			if (i == 0)
				x[i] = i * 1.5;
			w[i] = i * 1.5;
		}

		for (i = 0; i < n; i++) {
			if (i == 0) {
				z = (1.0 + alpha) * (3.0 + 0.92 * alpha) / (1.0 + 2.4 * n + 1.8 * alpha);
			}
			else if (i == 1) {
				z += (15.0 + 6.25 * alpha) / (1.0 + 0.9 * alpha + 2.5 * n);
			}
			else {
				ai = (double)i - 1;
				z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alpha /
					(1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * alpha);
			}
			for (its = 0; its < MAXIT; its++) {
				p1 = 1.0;
				p2 = 0.0;
				for (j = 0; j < n; j++) {
					p3 = p2;
					p2 = p1;
					p1 = ((2 * j + 1 + alpha - z) * p2 - (j + alpha) * p3) / ((double)j + 1);
				}
				pp = (n * p1 - (n + alpha) * p2) / z;
				z1 = z;
				z = z1 - p1 / pp;
				if (fabs(z - z1) <= EPS) break;
			}
			x[i] = z;
			w[i] = -exp(std::lgamma(alpha + n) - std::lgamma(double(n))) / (pp * n * p2);
			if (fabs(w[i]) <= MIN) break;
		}
	}

private:
		
};

