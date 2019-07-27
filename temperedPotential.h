#pragma once
#include <cmath>
#include "Gamma.h"   // macroscopic energy contribution

class temperedPotential
{
public:
	inline static double getTemperedRieszPotentialSphere(double alpha, double b,  double R0, double r) {
		
		if ( b == 0 ) return getRieszPotentialSphere( alpha, R0, r );

		if (alpha >= 2) {
			return 0;
		}

		if (r <= 0) {  // is normalized
			return -1.0;
		}
		
		Gamma gam;  // incomplete gamma function

		double pot;

		if ( r > R0 ) {
			pot = (b * (b * (-r * r + R0 * R0) * gam.gam(2. - alpha, b * (r - R0))
				+ b * (r - R0) * (r + R0) * gam.gam(2. - alpha, b * (r + R0))
				+ 2 * r * ( gam.gam(3. - alpha, b * (r - R0))
					- gam.gam(3. - alpha, b * (r + R0)) )) - gam.gam(4. - alpha, b * (r - R0)) +
				gam.gam(4. - alpha, b * (r + R0)) )	/ 
				(4. * b * r * (std::tgamma(3. - alpha) - gam.gam(3. - alpha, b * R0)) );
		}
		else {
			pot = (4. * b * r * std::tgamma(3. - alpha) + 
				b * b * (-r * r + R0 * R0) * gam.gam(2. - alpha, b * (-r + R0)) +
				b * (b * (r - R0) * (r + R0) * gam.gam(2. - alpha, b * (r + R0)) - 
					2. * r * (gam.gam(3. - alpha, b * (-r + R0)) + gam.gam(3. - alpha, b * (r + R0)))) 
				-gam.gam(4. - alpha, b * (-r + R0)) + gam.gam(4. - alpha, b * (r + R0))) /
				(4. * b * r * (std::tgamma(3. - alpha) - gam.gam(3. - alpha, b * R0)));
		}

		return -pot;  // normalized to V(r=0) = -1; 
			   		 	  	  
	}
	inline static double getRieszPotentialSphere(double alpha,  double R0, double r) {

		if ( alpha >= 3 || alpha <= 0) {
			return 0;
		}

		if (r <= 0) {  // is normalized
			return -1.0;
		}


		double a3 = 3.0 - alpha;
		double rm = r - a3 * R0;
		double rp = r + a3 * R0;

		double norm = 1/(r * 2 * (a3*a3 - 1) * pow(R0, 3.0 - alpha));

		double rm0 = pow(fabs(r - R0), a3);
		double rp0 = pow(r + R0, a3);
		
		return   (r > R0) ?  (rp0 * rm - rm0 * rp)*norm : (rp0 * rm + rm0 * rp)*norm;
	}

};

