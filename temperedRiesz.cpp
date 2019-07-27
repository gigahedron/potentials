// temperedRiesz.cpp : 
//

#include "pch.h"
#include <iostream>
#include <sstream>

#include <cmath>
#include <ctime>
#include <chrono>

#include "GHGaussL.h"   // Gauss Legendre abszissas and weights
#include "Gamma.h"   // macroscopic energy contribution
#include "temperedPotential.h"   // macroscopic energy contribution
#include "Laguerre.h"   // macroscopic energy contribution
#include "CompleteMath.h"   // macroscopic energy contribution


const int X = 9;
const int N = 10;


int main()
{
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	CompleteMath cm;
	
	cm.printGaussLegendre();
	cm.printLaguerre();
	
	
	double alpha = 1;
	double b = 0;
	double R0 = 10;
	double r = 10;
	double val = temperedPotential::getTemperedRieszPotentialSphere(alpha, b, R0, r);

	std::cout << "pot " << val << " qqq.";
	std::cout << std::endl;
	
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	return 0;
}

