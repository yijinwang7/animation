#include <string>
#include "Integrator.hpp"

class ModifiedMidpoint : public Integrator {
public:

    std::string getName() {
        return "modified midpoint";
    }

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 5, implmement the modified midpoint (2/3) method.
        // see also efficient memory management suggestion in provided code for the Midpoint method.
        //compute x(t) + 1/2*h*f(x)
        VectorXf& dpdt = VectorXf(n);
        derivs->derivs(t, p, dpdt);
        dpdt = p + dpdt * h * (2.0/3.0);

        //compute f(x(t) + 1/2*h*f(x))
        derivs->derivs(t + h * (2.0/3.0), dpdt, dpdt);
        //cpmpute x(t+h)
        pout = p + h * dpdt;
    }

};