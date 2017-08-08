/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// Random number generation stuff
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	mt19937 gen;

	// Set the number of particles
	num_particles = 100;

	// Create random particles around initial estimate
	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;
		particles.push_back(p);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Random number generation
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	mt19937 gen;
	
	// Move particles
	for (Particle& p : particles) {
		if (yaw_rate == 0) {
			p.x += velocity*delta_t*cos(p.theta);
			p.y += velocity*delta_t*sin(p.theta);
		} else {
			p.x += velocity / yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + dist_x(gen);
			p.y += velocity / yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + dist_y(gen);
			p.theta += yaw_rate*delta_t + dist_theta(gen);
		}
		// GPS uncertainty is used here because that's what's given...
		// Forum suggested angle normalization not needed
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// Some useful products to avoid repeat calculations
	double prob_function_denom = 1. / (2. * M_PI * std_landmark[0] * std_landmark[1]);
	double two_sigma_x_sq = 2. * std_landmark[0] * std_landmark[0];
	double two_sigma_y_sq = 2. * std_landmark[1] * std_landmark[1];

	// Clear weights
	weights.clear();

	for (Particle& p : particles) {
		// Reset the weight to 1
		p.weight = 1;

		for (LandmarkObs& obs : observations) {
			// Transform to global coordinates
			double x_m = p.x + cos(p.theta)*obs.x - sin(p.theta)*obs.y;
			double y_m = p.y + sin(p.theta)*obs.x + cos(p.theta)*obs.y;
			
			// Associate with the nearest map landmark
			double shortest_distance = 9999999999999999; // Initialize very high
			float x_lmk;
			float y_lmk;
			for (const auto& lmk : map_landmarks.landmark_list) {
				double d = dist(x_m, y_m, lmk.x_f, lmk.y_f);
				if (d < shortest_distance) {
					x_lmk = lmk.x_f;
					y_lmk = lmk.y_f;
					shortest_distance = d;
				}
			}

			// Calculate probablility of the observation, and multiply onto the weight
			double prob = prob_function_denom*exp(-((x_m - x_lmk)*(x_m - x_lmk) / two_sigma_x_sq + (y_m - y_lmk)*(y_m - y_lmk) / two_sigma_y_sq));
			p.weight *= prob;
		}

		// Populate weights vector - will help with resampling
		weights.push_back(p.weight);

		// Logging
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Generate the discrete distribution per particle weights
	discrete_distribution<> dist(weights.begin(), weights.end());
	mt19937 gen;

	// Resample
	vector<Particle> resampled_particles;
	for (int i = 0; i < num_particles; i++) {
		resampled_particles.push_back(particles[dist(gen)]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
