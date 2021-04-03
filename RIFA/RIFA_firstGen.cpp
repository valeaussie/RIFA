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
vector<double> F_founding_times(double n, double m0, double T0);
vector<double> F_detection_times(double n, double phi);
vector <double> F_vecsum(vector <double> a, vector <double> b);
double F_num_longjumps(double n, double r);
vector <double> F_direction(double n);
vector <double> F_short(double n, double r);
vector <double> F_long(double n, double r);
vector <double> F_x(double n, double parent, vector <double> d, vector <double> r);
vector <double> F_y(double n, double parent, vector <double> d, vector <double> r);

/* Simulate from the RIFA model and output the the observations. */
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
const float first_observation = 2000;
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
const double phi = 0.001;



int firstGen() {
	random_device rd;
	mt19937 generator(rd());


	//generation of the first set of nests from the first invader before the first set of observations arrive.
	//nnests are killed as soon as they are found
	//nests cannot produce new nests before they reach maturation (8 months).

	//time from first founding to first observation
	//during this times nests multiply undisturbed
	double T0 = first_observation - (first_nest_f + matur_time);
	//parameter of the Poisson distribution
	double lambda = rep_rate_daily * T0;


	vector <double> all_foundings;
	vector <double> all_detections;
	vector <double> all_x;
	vector <double> all_y;

	vector <double> temp_foundings;
	vector <double> temp_detections;
	vector <double> temp_x;
	vector <double> temp_y;

	do {
		//generate number of offsping for the first parent
		int num_off = F_num_offspring(lambda);

		//generate the time of funding of those offspring
		vector <double> temp_foundings = F_founding_times(num_off, rep_rate_daily, T0);

		//generate the time of detection of those offspring
		vector <double> detect = F_detection_times(num_off, phi);
		vector <double> temp_detections = F_vecsum(detect, temp_foundings);

		//generate the number of long jumps and short jumps
		int long_jumps = F_num_longjumps(num_off, long_jump_rate);
		int short_jumps = num_off - long_jumps;
		//generate radial direction
		vector <double> all_directions = F_direction(num_off);
		//generate short jumps
		vector <double> all_jumps = F_short(short_jumps, exp_jump_param);
		//generate long jumps
		vector <double> all_long_jumps = F_long(long_jumps, Levy_param);
		//vector <double> all_radial;
		all_jumps.insert(all_jumps.end(), all_long_jumps.begin(), all_long_jumps.end());
		//generate x and y coordinates
		vector <double> temp_x = F_x(num_off, first_nest_x, all_directions, all_jumps);
		vector <double> temp_y = F_y(num_off, first_nest_y, all_directions, all_jumps);

		cout << "test" << endl;

		all_foundings.insert(all_foundings.end(), temp_foundings.begin(), temp_foundings.end());
		all_detections.insert(all_detections.end(), temp_detections.begin(), temp_detections.end());
		all_x.insert(all_x.end(), temp_x.begin(), temp_x.end());
		all_y.insert(all_y.end(), temp_y.begin(), temp_y.end());


		//eliminate all observed nests that are not mature
		vector <double> foundings_ready;
		vector <double> detections_ready;
		vector <double> x_ready;
		vector <double> y_ready;
		for (int i = 0; i < temp_detections.size(); i++) {
			double detected = temp_detections[i];
			if (detected > first_observation) {
				double f_m = temp_foundings[i];
				double t_m = temp_detections[i];
				double x_m = temp_x[i];
				double y_m = temp_y[i];
				foundings_ready.push_back(f_m);
				detections_ready.push_back(t_m);
				x_ready.push_back(x_m);
				y_ready.push_back(y_m);
			}
		}
		vector <double> foundings_r_undetected;
		vector <double> detections_r_undetected;
		vector <double> x_r_undetected;
		vector <double> y_r_undetected;
		for (int i = 0; i < foundings_ready.size(); i++) {
			double mature = foundings_ready[i] + matur_time;
			if (mature < first_observation) {
				double f_d = foundings_ready[i];
				double t_d = detections_ready[i];
				double x_d = x_ready[i];
				double y_d = y_ready[i];
				foundings_r_undetected.push_back(f_d);
				detections_r_undetected.push_back(t_d);
				x_r_undetected.push_back(x_d);
				y_r_undetected.push_back(y_d);
			}
		}
		temp_foundings.clear();
		temp_foundings = foundings_r_undetected;
		temp_detections.clear();
		temp_detections = detections_r_undetected;
		temp_x.clear();
		temp_x = x_r_undetected;
		temp_y.clear();
		temp_y = y_r_undetected;

		cout << "test" << endl;
		cout << "temp_foundings " << endl;
		F_print_vector(temp_foundings);
	} while (temp_foundings.size() > 0);


	//add time of founding of first nest
	all_foundings.insert(all_foundings.begin(), first_nest_f);
	vector <double> found_detect = F_detection_times(1, phi);
	all_foundings.insert(all_foundings.end(), found_detect.begin(), found_detect.end());
	//add location of founding nest
	all_x.insert(all_x.begin(), first_nest_x);
	all_y.insert(all_y.begin(), first_nest_y);

	cout << "all_foundings " << endl;
	F_print_vector(all_foundings);


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
vector <double> F_founding_times(double n, double T0, double T1) {
	random_device rd;
	mt19937 generator(rd());
	vector <double> founding;
	uniform_real_distribution <double> unif(T0, T1);

	for (int i = 0; i < n; i++) {
		double f = unif(generator);
		founding.push_back(f);
	}
	return founding;
}

//F: simulates the detection times for n nests from a uniform distribution
vector <double> F_detection_times(double n, double phi) {
	random_device rd;
	mt19937 generator(rd());
	vector <double> detection;
	exponential_distribution <double> expo(phi);

	for (int i = 0; i < n; i++) {
		double f = expo(generator);
		detection.push_back(f);
	}
	return detection;
}
//F: summ the values of two vectors
vector <double> F_vecsum(vector <double> a, vector <double> b) {
	assert(a.size() == b.size());

	vector <double> sum;
	sum.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(), back_inserter(sum), plus<double>());

	return sum;
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
vector <double> F_direction(double n) {
	random_device rd;
	mt19937 generator(rd());
	vector < double > directions;
	uniform_real_distribution<double> unif_theta(0.0, 2 * M_PI);

	for (int i = 0; i < n; i++) {

		double f = unif_theta(generator);
		directions.push_back(f);
	}

	return directions;

}
double theta = 0;

//F: simulate the radial direction of the short jumps
vector <double> F_short(double n, double r) {
	random_device rd;
	mt19937 generator(rd());
	exponential_distribution<> exp_jump(r);
	vector < double > short_jumps;

	double short_rho = 0;

	for (int i = 0; i < n; i++) {
		double f = exp_jump(generator);
		short_jumps.push_back(f);
	}

	return short_jumps;

}

//F: simulate the radial direction of the long jumps
//there is no library for Levy distribution -> using inverse transform sampling
vector <double> F_long(double n, double r) {
	random_device rd;
	mt19937 generator(rd());

	double long_rho = 0;

	uniform_real_distribution<double> standardUnif(0.0, 1.0);

	vector < double > long_jumps;
	for (int i = 0; i < n; ++i) {
		double u = 1 - standardUnif(generator);
		double X = ((Levy_param) / ((1 / F_normalCDF(u)) * (1 / F_normalCDF(u))));
		long_jumps.push_back(X);
	}

	return long_jumps;
}

//F: find the x cartesian coordinates of each new nest
vector <double> F_x(double n, double parent, vector <double> d, vector <double> r) {
	vector < double > vec_x;
	for (int i = 0; i < n; ++i) {
		double x = 0;
		x = parent + d[i] * cos(r[i]);
		vec_x.push_back(x);
	}
	return vec_x;
}

//F: find the y cartesian coordinates of each new nest
vector <double> F_y(double n, double parent, vector <double> d, vector <double> r) {
	vector < double > vec_y;
	for (int i = 0; i < n; ++i) {
		double y = 0;
		y = parent + d[i] * sin(r[i]);
		vec_y.push_back(y);
	}
	return vec_y;
}


