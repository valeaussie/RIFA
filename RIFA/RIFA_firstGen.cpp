#define _USE_MATH_DEFINES
#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include <functional>
#include <assert.h>

using namespace std;

void F_print_vector(vector < double > value);
double F_normalCDF(double value);
double F_num_offspring(double m0);
double F_detection_times(double n, double phi);
double F_num_longjumps(double n, double r);
double F_direction(double n);
double F_short(double r);
double F_long(double r);
double F_x(double n, double parent, double d, double r);
double F_y(double n, double parent, double d, double r);
double F_lambda(double f, double t, double m, double T, double rep_rate_daily);
double F_founding_times(double n, double T, double t, double f, double fm);

/*Simulate from the RIFA model and output the the observations.*/
/*First invasion on Fisherman Island (Port of Brisbane)*/
//units of space in degrees and decimal minutes
//units of time in days
//latitude of the first nest
const double first_nest_x = -27.3541790;
//longitude of the first nest
const double first_nest_y = 153.1930033;
//time of founding of the first nest
const float first_nest_f = 0;
//time of first observation (in days)
const float observ_time = 2000;
//reproductive rate number of nests per month
const double rep_rate = 0.25;
//reproductive rate number of nests per day
const double rep_rate_daily = rep_rate * 12 / 365;
//probability of a long jump (this value is a guess at the moment)
const double long_jump_rate = 0.1;
//parameter for the Levy distribution (this value is a guess at the moment)
const double Levy_param = 100;
//parameter for the exponential distribution for the short distance jumps (this value is a guess at the moment)
const double exp_jump_param = 1;
//maturation time for a nest in days (average 30.42 days in a month)
const double matur_time = 8 * 30.42;
//parameter of the exponential distribution for the times of detection
const double phi = 0.00001;

//nests are killed as soon as they are found
//nests cannot produce new nests before they reach maturation (8 months).

int firstGen() {
	random_device rd;
	mt19937 generator(rd());

	//this is the catalogue of all nests
	vector <double> all_foundings = { first_nest_f };
	vector <double> all_detections = { 20000 };
	vector <double> all_x = { -27.3541790 };
	vector <double> all_y = { 153.1930033 };

	//this is the G catalogue
	vector <double> t_foundings = all_foundings;
	vector <double> t_detections = all_detections;
	vector <double> t_x = all_x;
	vector <double> t_y = all_y;

	//this is the O catalogue
	vector <double> t_t_foundings;
	vector <double> t_t_detections;
	vector <double> t_t_x;
	vector <double> t_t_y;

	bool notempty = true;

	while (notempty) {
		for (int i = 0; i < t_foundings.size(); i++) {

			//calculate the parameter of the Poisson distribution
			double lambda = -1;
			lambda = F_lambda(t_foundings[i], t_detections[i], matur_time, observ_time, rep_rate_daily);
			if (lambda > 1) {

				//generate number of offsping for the parent i
				int num_off = F_num_offspring(lambda);

				//generate the number of long jumps and short jumps
				int n_long_jumps = F_num_longjumps(num_off, long_jump_rate);
				int n_short_jumps = num_off - n_long_jumps;


				//create a vector for the radial directions
				vector <double> radial;
				for (int j = 0; j < n_short_jumps; j++) {
					//generate short jumps
					double short_jumps = F_short(exp_jump_param);
					radial.push_back(short_jumps);
				}
				for (int j = 0; j < n_long_jumps; j++) {
					//generate long jump
					double long_jumps = F_long(Levy_param);
					radial.push_back(long_jumps);
				}
				
				for (int j = 0; j < num_off; j++) {
					//generate the time of funding of the offspring
					double f = F_founding_times(num_off, observ_time, t_detections[i], t_foundings[i], matur_time);
					t_t_foundings.push_back(f);
					all_foundings.push_back(f);

					//generate the time of detection of the offspring
					double t = F_detection_times(num_off, phi);
					t_t_detections.push_back(t);
					all_detections.push_back(t);

					//generate radial direction
					double direction = F_direction(num_off);

					//generate x and y coordinates
					double x = F_x(num_off, t_x[i], direction, radial[j]);
					t_t_x.push_back(x);
					all_x.push_back(x);

					double y = F_y(num_off, t_y[i], direction, radial[j]);
					t_t_y.push_back(y);
					all_y.push_back(y);
				}
			}
		}

		t_foundings.clear();
		t_detections.clear();
		t_x.clear();
		t_y.clear();

		t_foundings = t_t_foundings;
		t_detections = t_t_detections;
		t_x = t_t_x;
		t_y = t_t_y;

		t_t_foundings.clear();
		t_t_detections.clear();
		t_t_x.clear();
		t_t_y.clear();

		if (t_foundings.size() == 0)  notempty = false; 
	}	

	//Output CSV
	ofstream outFile("./NestsSym.dat");
	outFile << "Foundings, Detections, x, y" << endl;
	for (size_t lin = 0; lin << t_t_foundings.size(); lin++) {
		outFile << "t_t_foundings[lin], t_t_detections[lin], t_t_x[lin], t_t_y[lin]" << endl;
	}
	outFile.close();

	return 0;
}



//*****************************************************
//*													  *
//*				 FUNCTIONS DEFINITION				  *
//*													  *
//*****************************************************


//F: prints a vector of doubles
void F_print_vector(vector < double > value) {
	for (const double x : value) cout << x << ' ';
	cout << endl;
}

//F: calculates the CDF of a standard normal distribution
double F_normalCDF(double value) {
	return 0.5 * erfc(-value * M_SQRT1_2);
}

//F: generate the time T0 from mature to first observation
//the calculate parameters of the Poisson distributions
double F_lambda(double f, double t, double m, double T, double rep_rate_daily) {
	double T0 = 0;
	double lambda = -1;

	if (f < t && T > (f + m)) {
		if (T > t) {
			T0 = t - (f + m);
		}
		else {
			T0 = T - (f + m);
		}
		if (T0 > 0) {
			lambda = rep_rate_daily * T0;
		}
	}
	return lambda;
}

//F: simulates the number of nests from a poisson distribution
double F_num_offspring(double m0) {
	random_device rd;
	mt19937 generator(rd());
	poisson_distribution<int> pois(m0);
	int num_offsp = pois(generator);
	if (num_offsp == 0) {
		num_offsp = 1;
	}
	return num_offsp;
}

//F: simulates the founding times for n nests from a uniform distribution
double F_founding_times(double n, double T, double t, double f, double fm) {
	random_device rd;
	mt19937 generator(rd());
	double founding;
	double T0 = f + fm;
	double T1 = 0;
	if (t > T) {
		T1 = T;
	}
	else { T1 = t; }
	uniform_real_distribution <double> unif(T0, T1);
	founding = unif(generator);
	return founding;
}

//F: simulates the detection times for n nests from a uniform distribution
double F_detection_times(double n, double phi) {
	random_device rd;
	mt19937 generator(rd());
	exponential_distribution <double> expo(phi);
	double detection = expo(generator);
	return detection;
}


//F: simulates number of long distance jumps
double F_num_longjumps(double n, double r) {
	random_device rd;
	mt19937 generator(rd());
	binomial_distribution <> binom(n, r);
	double num_long_jump = binom(generator);
	return num_long_jump;
}

//F: simulate the angular direction of the short jumps (short and long)
double F_direction(double n) {
	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> unif_theta(0.0, 2 * M_PI);
	double directions = unif_theta(generator);
	return directions;
}

//F: simulate the radial direction of the short jumps
double F_short(double r) {
	random_device rd;
	mt19937 generator(rd());
	exponential_distribution<> exp_jump(r);
	double short_jumps = exp_jump(generator);
	return short_jumps;
}

//F: simulate the radial direction of the long jumps
//there is no library for Levy distribution -> using inverse transform sampling
double F_long(double r) {
	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> standardUnif(0.0, 1.0);
	double u = 1 - standardUnif(generator);
	double long_jumps = ((Levy_param) / ((1 / F_normalCDF(u)) * (1 / F_normalCDF(u))));
	return long_jumps;
}

//F: find the x cartesian coordinates of each new nest
double F_x(double n, double parent, double d, double r) {
	double x = 0;
	x = parent + d * cos(r);
	return x;
}

//F: find the y cartesian coordinates of each new nest
double F_y(double n, double parent, double d, double r) {
	double y = 0;
	y = parent + d * sin(r);
	return y;
}

//F: this function writes a dat file 




