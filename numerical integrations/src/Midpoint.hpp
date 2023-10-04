#include <string>
#include "Integrator.hpp"

class Midpoint : public Integrator {
public:
    std::string getName() {
        return "midpoint";
    }

    float* tmp;
    int tmplength;

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 4, implement midpoint method

        // You will probably want a temporary array in this method 

        //compute x(t) + 1/2*h*f(x)
        VectorXf& dpdt = VectorXf(n);
        derivs->derivs(t, p, dpdt);
        dpdt = p + dpdt * h * 0.5;

        //compute f(x(t) + 1/2*h*f(x))
        derivs->derivs(t+h*0/5, dpdt, dpdt);
        //cpmpute x(t+h)
        pout = p + h * dpdt;

    }

};
