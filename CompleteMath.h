#pragma once

#include "pch.h"
#include <iostream>
#include <sstream>

#include <cmath>
#include <ctime>
#include <chrono>

#include "GHGaussL.h"			// Gauss Laguerre abszissas and weights
#include "Gamma.h"				// macroscopic energy contribution
#include "temperedPotential.h"  // macroscopic energy contribution
#include "Laguerre.h"			// macroscopic energy contribution

class CompleteMath
{

private:
	const static long N = 10;
	const static long X = 18;  // nGauss


	// Gauss Laguerre abscissas and weights
	
	long nGauss = X;
	double alphaGauss = 1.0;
	std::vector<double> abscissaArray;
	std::vector<double> weightArray;

	// Legendre

	long nLaguerre = N;
	long jLaguerre = 0;
	std::vector<double> LaguerreN;
	std::vector<double> LaguerreNorm;
	double matLaguerreL[N][X];



public:
	CompleteMath() {
		setGaussLegendre(nGauss, 1.0); // start with n=18, alpha = 1;
		setLaguerreArray(nLaguerre, 0);
	};
	~CompleteMath() {};

	inline double setGaussLegendre (long n, double alpha) {
		
		// initialization of abscissas and weights
		nGauss = n;
		alphaGauss = alpha;
		abscissaArray.empty();
		weightArray.empty();

		abscissaArray.reserve(n);
		weightArray.reserve(n);

		GHGaussL::getList(n, alpha, abscissaArray, weightArray);  // generate abscissas and weights
		
		return nGauss;

	}

	inline void printGaussLegendre() {
		
		
		std::cout << "nGauss " << nGauss << std::endl;
		std::cout << "alphaGauss " << alphaGauss << std::endl;
	
		for (long i = 0; i < nGauss; i++)
			std::cout << "i " << i << " : " << std::setprecision(18) << abscissaArray[i] << " " << weightArray[i] << std::endl;
	
	}

	inline void setLaguerreArray(long N, long J, bool normalize = true) {
	
		// initialization 
		nLaguerre = N;
		jLaguerre = J;

		LaguerreN.empty();
		LaguerreNorm.empty(); 

		LaguerreN.reserve(nLaguerre);
		LaguerreNorm.reserve(nLaguerre);

		Laguerre::getNorm(nLaguerre, jLaguerre + 0.5, LaguerreNorm);  // for a given angular momentum J 

		for (long j = 0; j < nGauss; j++) {
			double x = abscissaArray[j];
			Laguerre::generateList(nLaguerre, jLaguerre + 0.5, x, LaguerreN, normalize, LaguerreNorm);
			for (long i = 0; i < nLaguerre; i++) {
				matLaguerreL[i][j] = LaguerreN[i];
			}
		}

	}

	inline void printLaguerre() {
		std::cout << "nLaguerre " << nLaguerre << std::endl;
		std::cout << "X " << nGauss << std::endl;

		for (long i = 0; i < nLaguerre; i++) {
			for (long j = 0; j < nGauss; j++) {
				std::cout << "Lag_ " << i << "("<< abscissaArray[j] << ")= " << matLaguerreL[i][j] << std::endl;
			}
		}
	}


};

