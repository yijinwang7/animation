#include <string>
#include "Integrator.hpp"

class RK4 : public Integrator {
public:
    std::string getName() {
        return "RK4";
    }

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 6, implement the RK4 integration method
        // see also efficient memory management suggestion in provided code for the Midpoint method.
        VectorXf& k1 = VectorXf(n);
        VectorXf& k2 = VectorXf(n);
        VectorXf& k3 = VectorXf(n);
        VectorXf& k4 = VectorXf(n);
        //k1
        derivs->derivs(t, p, k1);
        k1 = k1 * h;
        //k2
        k2 = p + k1 * 0.5;
        derivs->derivs(t + h * 0.5, k2, k2);
        k2 = k2 * h;
        //k3
        k3 = p + k2 * 0.5;
        derivs->derivs(t + h * 0.5, k3, k3);
        k3 = k3 * h;
        //k4
        k4 = p + k3;
        derivs->derivs(t + h, k4, k4);
        k4 = k4 * h;

        pout = p + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
};