#pragma once
#include "pch.h"
#include <algorithm>  

struct Gauleg18 {
	static const long ngau = 18;
	static const double y[18];
	static const double w[18];
};
const double Gauleg18::y[18] = { 0.0021695375159141994,
0.011413521097787704,0.027972308950302116,0.051727015600492421,
0.082502225484340941, 0.12007019910960293,0.16415283300752470,
0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
0.87126389619061517, 0.95698180152629142 };
const double Gauleg18::w[18] = { 0.0055657196642445571,
0.012915947284065419,0.020181515297735382,0.027298621498568734,
0.034213810770299537,0.040875750923643261,0.047235083490265582,
0.053244713977759692,0.058860144245324798,0.064039797355015485,
0.068745323835736408,0.072941885005653087,0.076598410645870640,
0.079687828912071670,0.082187266704339706,0.084078218979661945,
0.085346685739338721,0.085983275670394821 };
struct Gamma : Gauleg18 {
	static const long ASWITCH = 100;
	static  constexpr double EPS = 1.E-32;
	static  constexpr double FPMIN = 1.E-32;
	double gln;

	double gammp(const double a, const double x) {   // P(a,x)
		if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
		if (x == 0.0) return 0.0;
		else if ((long)a >= ASWITCH) return gammpapprox(a, x, 1);
		else if (x < a + 1.0) return gser(a, x);
		else return 1.0 - gcf(a, x);
	}

	double gammq(const double a, const double x) {
		if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
		if (x == 0.0) return 1.0;
		else if ((long)a >= ASWITCH) return gammpapprox(a, x, 0);
		else if (x < a + 1.0) return 1.0 - gser(a, x);
		else return gcf(a, x);
	}

	double gam(const double a, const double x) {  // mathematica definition
		return std::tgamma(a)*gammq(a, x);
	}



	double gser(const double a, const double x) {
		double sum, del, ap;
		gln = lgamma(a);
		ap = a;
		del = sum = 1.0 / a;
		for (;;) {
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS) {
				return sum * exp(-x + a * log(x) - gln);
			}
		}
	}

	double gcf(const double a, const double x) {
		long i;
		double an, b, c, d, del, h;
		gln = lgamma(a);
		b = x + 1.0 - a;
		c = 1.0 / FPMIN;
		d = 1.0 / b;
		h = d;
		for (i = 1;; i++) {
			an = -i * (i - a);
			b += 2.0;
			d = an * d + b;
			if (fabs(d) < FPMIN) d = FPMIN;
			c = b + an / c;
			if (fabs(c) < FPMIN) c = FPMIN;
			d = 1.0 / d;
			del = d * c;
			h *= del;
			if (fabs(del - 1.0) <= EPS) break;
		}
		return exp(-x + a * log(x) - gln) * h;
	}

	double gammpapprox(double a, double x, long psig) {
		long j;
		double xu, t, sum, ans;
		double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
		gln = lgamma(a);
		if (x > a1) xu = std::max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
		else xu = std::max(0., std::min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
		sum = 0;
		for (j = 0; j < ngau; j++) {
			t = x + (xu - x) * y[j];
			sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
		}
		ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
		return (psig ? (x > a1 ? 1.0 - ans : -ans) : (x > a1 ? ans : 1.0 + ans));
	}

};

