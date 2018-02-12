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



#define THRESHOLD 0.001

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//set the number of particles
	num_particles = 100;

	default_random_engine gen;

	//create guassian distribution
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	//initialize values to all the particles
	for(int i = 0; i < num_particles; i++){
		Particle P;
		P.id = i;
		P.x = dist_x(gen);
		P.y = dist_y(gen);
		P.theta = dist_theta(gen);
		P.weight = 1.0;
		particles.push_back(P);
	}

	is_initialized = true;
	cout << "completed Init" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	double new_theta;
	double new_x;
	double new_y;

	for(int i = 0; i < num_particles; i++){

		// check if the yaw is not zero
		if(fabs(yaw_rate) > THRESHOLD){
			//update measurements to particles
			new_theta = particles[i].theta + yaw_rate * delta_t;
			new_x = particles[i].x + (velocity / yaw_rate) * (sin(new_theta) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(new_theta));
		}
		else{
			new_theta = particles[i].theta;
			new_x = particles[i].x + velocity * delta_t * cos(new_theta);
			new_y = particles[i].y + velocity * delta_t * sin(new_theta);
		}

		//noise generated to be added
		normal_distribution<double> dist_x(new_x, std_pos[0]);
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		//adding noise
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
//cout << "completed prediction" << endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
/* Not using
	unsigned int num_observations = observations.size();
	unsigned int num_predictions = predicted.size();

	for(int i = 0; i < num_observations; i++){
		double min_distance = numeric_limits<double>::max();
		int map_id = -1;

		for(int j = 0; j < num_predictions; j++){

			double x_distance = observations[j].x - predicted[j].x;
			double y_distance = observations[j].y - predicted[j].y;

			double distance = x_distance * x_distance + y_distance * y_distance;

			if(distance < min_distance){
				min_distance = distance;
				map_id = predicted[j].id;
			}
		}
		observations[i].id = map_id;
	}
	//cout << "completed data association" << endl;

*/
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	//cout << "number of part: " << num_particles << endl;
	//for(int i = 0; i < num_particles; i++){
	//			cout << particles[i].weight << " ";
	//		}
	//		cout << endl << endl;

	LandmarkObs mappedObservation;
	LandmarkObs inRange;

	weights.clear();

	for(int i = 0; i < num_particles; i++){

		double part_weight = 1.0;
		double part_x = particles[i].x;
		double part_y = particles[i].y;
		double part_theta = particles[i].theta;

		for(int j = 0; j < observations.size(); j++){

			double mapped_x = cos(part_theta) * observations[j].x - sin(part_theta) * observations[j].y + part_x;
			double mapped_y = sin(part_theta) * observations[j].x + cos(part_theta) * observations[j].y + part_y;
			mappedObservation = LandmarkObs { observations[j].id, mapped_x, mapped_y };

			for (int k = 0; k < int(map_landmarks.landmark_list.size()); k++) {
				double min_dist;
				double distance = dist(mappedObservation.x, mappedObservation.y,
						map_landmarks.landmark_list[k].x_f,
						map_landmarks.landmark_list[k].y_f);

				if (k == 0) {
					min_dist = distance;
					inRange.id = map_landmarks.landmark_list[k].id_i;
					inRange.x = map_landmarks.landmark_list[k].x_f;
					inRange.y = map_landmarks.landmark_list[k].y_f;
				} else if (distance < min_dist) {
					min_dist = distance;
					inRange.id = map_landmarks.landmark_list[k].id_i;
					inRange.x = map_landmarks.landmark_list[k].x_f;
					inRange.y = map_landmarks.landmark_list[k].y_f;
				}
			}

			double diff_x = mappedObservation.x - inRange.x;
			double diff_y = mappedObservation.y - inRange.y;

			double gauss_norm = 2 * M_PI * std_landmark[0] * std_landmark[1];
			double a = 2 * std_landmark[0] * std_landmark[0];
			double b = 2 * std_landmark[1] * std_landmark[1] ;
			double exponent = (diff_x * diff_x / a) + (diff_y * diff_y / b);
			double observed_weight = exp( -(exponent) ) / gauss_norm;

			part_weight *= observed_weight;
		}

		particles[i].weight = part_weight;

		weights.push_back(part_weight);
	}

	//cout << endl;
	//for(int i = 0; i < num_particles; i++){
	//		cout << particles[i].weight << " ";
	//}
	//	cout << endl;

	//cout << "completed weight update" << endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


	  std::vector<Particle> resampled_particle;
	  int index;

	  // Initializes discrete distribution function
	  std::random_device rd;
	  std::mt19937 gen(rd());
	  std::discrete_distribution<int> weight_distribution(weights.begin(), weights.end());

	  for (int i = 0; i < num_particles; i++) {
	    index = weight_distribution(gen);
	    resampled_particle.push_back(particles[index]);
	  }
	  particles = resampled_particle;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
