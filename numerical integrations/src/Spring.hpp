#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <Eigen/Dense>
using Eigen::MatrixXf;
using Eigen::VectorXf;

#include "Particle.hpp"

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
class Spring {

public:

    Particle* p1;
    Particle* p2;

    /** Spring stiffness, sometimes written k_s in equations */
    float k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    float c = 1;
    /** Rest length of this spring */
    double l0 = 0;

    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    Spring(Particle* p1, Particle* p2) {
        this->p1 = p1;
        this->p2 = p2;
        recomputeRestLength();
        p1->springs.push_back(this);
        p2->springs.push_back(this);
    }

    /**
     * Computes and sets the rest length based on the original position of the two particles
     */
    void recomputeRestLength() {
        glm::vec2 diff = p1->p0 - p2->p0;
        l0 = sqrt(diff.x * diff.x + diff.y * diff.y);
    }

    /**
     * Applies the spring force by adding a force to each particle
     */
    void apply() {
        // TODO: Objective 1, FINISH THIS CODE!
        
            //Set the line vector
            glm::vec2 line = this->p2->p - this->p1->p;
            //get the distance of the two particles
            float d_line = this->p1->distance(this->p2->p.x, this->p2->p.y);
            //Set up the line velocity
            glm::vec2 velocity = this->p2->v - this->p1->v;
            // scaling of line
            //float s = (k * (d_line  - l0) + c * glm::dot(velocity, line) / (l0 * d_line)) / (d_line);
            
            float s = k * (d_line-l0);
            s = s - c * glm::dot(velocity,line/d_line);
            //add forces
            if (!p1->pinned) {
                this->p1->addForce(s * line/d_line);
            }
            if (!p2->pinned) {
                this->p2->addForce(-s * line/d_line);
            }

            
        
        
        
        
    }

    /** TODO: the functions below are for the backwards Euler solver */

    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * (This function is something you might use for a backward Euler integrator)
     * @param f
     */
    void addForce(VectorXf& f) {
        // TODO: Objective 8, FINISH THIS CODE for backward Euler method (probably very simlar to what you did above)

    }

    /**
     * Adds this springs contribution to the stiffness matrix
     * @param dfdx
     */
    void addDfdx(MatrixXf& dfdx) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration

    }

    /**
     * Adds this springs damping contribution to the implicit damping matrix
     * @param dfdv
     */
    void addDfdv(MatrixXf& dfdv) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward Euler integration

    }

};