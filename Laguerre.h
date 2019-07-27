#pragma once
#include "pch.h"
#include <vector>
#include <cmath> 

struct Laguerre  {
	
	inline static double generateList(const long n, const double alpha, const double x, 
		std::vector<double>& LaguerreN, 
		bool normalize = false, const std::vector<double>& vecNorm = std::vector<double>()) {   // P(a,x)
		std::cout << "nnn " << n << std::endl;
		LaguerreN[0] = 1.0 ;
		LaguerreN[1] = (1.0 -x + alpha);
	
		for (long i = 1; i < n-1; i++) {
			LaguerreN[i + 1] = ((2.0 * i + 1.0 + alpha - x) * LaguerreN[i] - (i + alpha) * LaguerreN[i - 1]) / (i + 1.0);
		}
		if (normalize) {
			for (long i = 0; i < n; i++) {
				std::cout << "unnormalized Lag_" << i << " " << LaguerreN[i] << std::endl;
				LaguerreN[i] *= vecNorm[i];
				std::cout << "normalized Lag_" << i << " " << LaguerreN[i] << std::endl;
			}
		}
		return x;
	}
	
	inline static double getNorm(const long n, const double alpha, std::vector<double>& vecNorm) {   // P(a,x)
		std::cout << "norm " << n << std::endl;
		double sqrt2 = 1.4142135623730950488016887242097; // sqrt(2)
		double J = alpha - 0.5;
		
		for (long i = 0; i < n; i++) {
			vecNorm[i] = sqrt2 * exp( (lgamma(1.0 + i) - lgamma(1.5 + J + i) )* 0.5);
			std::cout << "norm  Lag_" << i << " " << vecNorm[i] << std::endl;
		}
		return vecNorm[n-1];
	}

	

};

