#pragma once
#include "Particle.hpp"
#include "Spring.hpp"
#include <algorithm>

class RobustCCD {
public:
	
	bool collision;
	bool useRepulsion;
	int maxIterations;
	double restitutionValue;
	double thresholdDistance;
	double eps;
	double repulsionStiffness;

	/** TODO: number of iterations in the last CCD processing loop, to keep an eye on how tricky the problem is */
	int iterations;
	/** TODO: keep track of highest count needed so far */
	int highestIterations = 0;

	RobustCCD::RobustCCD() :
		eps(0.000001),
		//eps(0.001),
		collision(true),
		useRepulsion(true),
		//maxIterations(100),
		//restitutionValue(0.0),
		maxIterations(1),
		restitutionValue(1.0),
		thresholdDistance(2),
		iterations(0),
		repulsionStiffness(100) {}

	/**
	 * Try to deal with contacts before they happen
	 * @param h
	 * @param system
	 */
	void applyRepulsion(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs) {
		//TODO: apply repulsion on all particles

		// - check for point-edge proximity between all particle-edge pairs
		// - find the normal
		// - take care to deal with segment end points carefully
		// - compute an appropriate  impulse if the distance is less than H
		// - make sure your impulse is in the correct direction!
		// - don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
		// - distribute impulse to the three particles involved in the appropriate manner
		for (Particle* p0 : *particles) {

			for (Spring* s : *springs) {
				if (p0 == s->p1 || p0 == s->p2) continue;
				glm::vec2 diff2 = p0->p - s->p1->p;
				glm::vec2 diff1 = s->p2->p - s->p1->p;
				glm::vec2 norm = normalize(diff1);
				double x2 = -norm.y;
				double y2 = norm.x;
				norm.x = x2;
				norm.y = y2;
				//if (glm::dot(diff1, diff2) > 0) norm = -norm;
				double distance = glm::dot(diff2, norm);
				if (distance  < thresholdDistance && distance >0) {
					double alpha = (diff1.x * diff2.x + diff1.y * diff2.y) / (diff1.x * diff1.x + diff1.y * diff1.y);
					if(alpha > -eps && alpha < 1 + eps){
						double p0Mass = p0->mass;
						double p1Mass = s->p1->mass;
						double p2Mass = s->p2->mass;
						//if (p0->pinned) p0Mass = INFINITY;
						//if (s->p1->pinned) p1Mass = INFINITY;
						//if (s->p2->pinned) p2Mass = INFINITY;
						//first impulse
						glm::vec2 rel_vol = alpha * s->p1->v + (1 - alpha) * s->p2->v - p0->v;
						double normalComponentRel_vel = glm::dot(norm, rel_vol);
						double impulse1 = 0.25*normalComponentRel_vel / 2;
						p0->v.x = p0->v.x + impulse1 * norm.x;
						s->p1->v.x = s->p1->v.x + (-alpha)*impulse1 * norm.x;
						s->p2->v.x = s->p2->v.x + (alpha-1)*impulse1 * norm.x;
						p0->v.y = p0->v.y + impulse1 * norm.y;
						s->p1->v.y = s->p1->v.y + (-alpha)*impulse1 * norm.y;
						s->p2->v.y = s->p2->v.y + (alpha-1)*impulse1 * norm.y;
						
						//second impulse

						glm::vec2 rel_pos = p0->p - alpha * s->p1->p - (1 - alpha) * s->p2->p;
						double d = thresholdDistance - glm::dot(rel_pos, norm);
						if (normalComponentRel_vel <  0.1 * d / h) {
							double k = h * repulsionStiffness * d;
							double impulseP0 = (0.1 * d / h) - normalComponentRel_vel;
							double impulseP1 = (0.1 * d / h) - normalComponentRel_vel;
							double impulseP2 = (0.1 * d / h) - normalComponentRel_vel;
							if (k < impulseP0 * p0Mass) impulseP0 = k / p0Mass;
							if (k < impulseP1 * p1Mass) impulseP1 = k / p1Mass;
							if (k < impulseP2 * p2Mass) impulseP2 = k / p2Mass;
							/*p0->v = p0->v + impulseP0*norm;
							s->p1->v = s->p1->v + impulseP1*norm;
							s->p2->v = s->p2->v + impulseP2*norm;*/
							p0->v.x -= impulseP0 * norm.x;
							s->p1->v.x -= -alpha*impulseP1 * norm.x;
							s->p2->v.x -= (alpha-1)*impulseP2 * norm.x;
							p0->v.y -= impulseP0 * norm.y;
							s->p1->v.y -= -alpha*impulseP1 * norm.y;
							s->p2->v.y -= (alpha-1)*impulseP2 * norm.y;
						}
						useRepulsion = true;
					}
					
					
					
				}
				
				
			}
		}
		
	}

	/**
	* Checks all collisions in interval t to t+h
	* @param h
	* @param system
	* @return true if all collisions resolved
	*/
	bool checkCollision(double h, std::vector<Particle*>* particles, std::vector<Spring*>* springs) {
		// TODO: Most of your assignment code will go here
		// - put your iterative solving code here!
		// - you can use the nextCollision function below to deal with an individual particle-edge pair

		// For each particle-edge pair, find the roots for when the three particles are
		// co-linear, and then pick the first root on (0,h] which corresponds to an actually collision.
		// compute a collision response.  That is, compute an appropriate collision normal, compute the 
		// impulse, and then apply the impulse to the associated particles.  Be sure to deal with pinning 
		// constraints!


		iterations = 0;
		
		while (iterations < maxIterations) {
			collision = false;

			for (Particle* p : *particles) {

				for (Spring* s : *springs) {

					if (p == s->p1 || p == s->p2) continue;
					 
					if (nextCollision(h, restitutionValue, p, s->p1, s->p2)) collision = true;
				}
			}

			if (!collision) {
				break;
			}
			
			iterations += 1;
			
		}
		return true;


		
	}

	/**
	* Processes next collision on (0,h] if it exists, between p0 and segment connecting p1 and p2.
	*
	* Watch out for
	* - normal direction
	* - epsilon tests
	* - segment connecting p1 and p2 at zero length! (unlikely?)
	*
	* @param h				timestep
	* @param restitution	bouncyness parameter
	* @param p0				particle
	* @param p1				line segment end point
	* @param p2				line segment end point
	* @return true if collision processed
	*/
	bool RobustCCD::nextCollision(double h, double restitution, Particle* p0, Particle* p1, Particle* p2) {
		//TODO: 
		// * - finds roots
		// * - checks that p0 actually falls on segment
		// * - processes collision by computing and applying impulses
		// * - returns true if collision processed
		
		//form the equation a*t^2+b*t+c=0
		double Adot1 = p1->v.x;
		double Adot2 = p1->v.y;
		double Bdot1 = p2->v.x;
		double Bdot2 = p2->v.y;
		double A1 = p1->p.x;
		double A2 = p1->p.y;
		double B1 = p2->p.x;
		double B2 = p2->p.y;
		double C1 = p0->p.x;
		double C2 = p0->p.y;
		double Cdot1 = p0->v.x;
		double Cdot2 = p0->v.y;
		double a = Adot1 * Bdot2 - Adot2 * Bdot1 - Adot1 * Cdot2 + Adot2 * Cdot1 + Bdot1 * Cdot2 - Bdot2 * Cdot1;
		double b = A1 * Bdot2 - A2 * Bdot1 + Adot1 * B2 - Adot2 * B1 - A1 * Cdot2 + A2 * Cdot1 - Adot1 * C2 + Adot2 * C1 +B1 * Cdot2 - B2 * Cdot1 + Bdot1 * C2 - Bdot2 * C1;
		double c = A1 * B2 - A2 * B1 - A1 * C2 + A2 * C1 + B1 * C2 - B2 * C1;
		std::vector<double> root;
		if (a == 0) {
			if (b != 0) {
				root.push_back(-c / b);
			}

			else {
			    root.push_back(0.0); // treat this like area = 1/2 * c and check the boundary < eps
			}
		}
		else {
			double delta = b * b - 4 * a * c;
			if (delta == 0) {
				root.push_back((-b) / (2 * a));
			}
			else if (delta > 0) {
				root.push_back((-b - sqrt(delta)) / (2 * a));
				root.push_back((-b + sqrt(delta)) / (2 * a));
			}
			
		}

		for (double t : root) {
			if (t > -eps && t <= h + eps ) {
				glm::vec2 newP0 = p0->p + p0->v * t;
				glm::vec2 newP1 = p1->p + p1->v *t;
				glm::vec2 newP2 = p2->p +p2->v *t;
				glm::vec2 diff2 = newP0-newP1;
				glm::vec2 diff1 = newP2-newP1;
				double boundary = (a * t * t + b * t + c) / sqrt(glm::dot(diff1, diff1));
				double alpha = (diff1.x * diff2.x + diff1.y * diff2.y) / (diff1.x * diff1.x + diff1.y * diff1.y);
				if (alpha >= -eps && alpha <= 1 + eps && abs(boundary) < eps) {
					
					/*p1->pinned = true;
					p2->pinned = true;
					p0->pinned = true;*/
					
					
					//compute unit normal
					glm::vec2 norm = normalize(newP2-newP1 );
					double x2 = -norm.y;
					double y2 = norm.x;
					norm.x = x2;
					norm.y = y2;
					//compute the impulse

					// pinned particles have infinite mass

					double p0Mass = p0->mass;
					double p1Mass = p1->mass;
					double p2Mass = p2->mass;

					if (p0->pinned) {
						p0Mass = INFINITY;
					}
				
					if (p1->pinned) {
						p1Mass = INFINITY;
					}
				
					if (p2->pinned) {
						p2Mass = INFINITY;
					}
				
					glm::vec2 rel_vol = alpha * p1->v + (1 - alpha) * p2->v - p0->v;
					double pp1 = alpha * alpha / p1Mass;
					double pp2 = (1 - alpha) * (1 - alpha) / p2Mass;
					double impulse = (1 + restitutionValue) * glm::dot(rel_vol, norm) / (pp1 + pp2 + (1 / p0Mass));
					p0->v.x = p0->v.x + impulse * norm.x / p0Mass;
					p0->v.y = p0->v.y + impulse * norm.y / p0Mass;

					p1->v.x = p1->v.x + alpha*(-impulse * norm.x / p0Mass);
					p1->v.y = p1->v.y + alpha*(-impulse * norm.y / p0Mass);

					p2->v.x = p2->v.x + (alpha-1)*impulse * norm.x / p0Mass;
					p2->v.y = p2->v.y + (alpha-1)*impulse * norm.y / p0Mass;

					//p0->v = p0->v + impulse * norm / p0Mass;
					//p1->v = p1->v + alpha * (-impulse * norm / p1Mass);
					//p2->v = p2->v + (1-alpha) * (-impulse * norm / p2Mass);
					return true;
				}
			}
		}
		return false;
		
	}


};
