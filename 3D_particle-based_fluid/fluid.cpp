//3d version
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <vector>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

// "Particle-Based Fluid Simulation for Interactive Applications" by Müller et al.
// solver parameters
const static float M_PI = 3.1415926535f;
const static Vector3d G(0.f, -10.f, 0.f);   // external (gravitational) forces
const static float REST_DENS = 300.f;  // rest density
const static float GAS_CONST = 2000.f; // const for equation of state
const static float H = 16.f;		   // kernel radius
const static float HSQ = H * H;		   // radius^2 for optimization
const static float MASS = 2.5f;		   // assume all particles have the same mass
const static float VISC = 200.f;	   // viscosity constant
const static float DT = 0.0007f;	   // integration timestep
const float sigma = 0.0728f;           // surface tension constant
const float epsilon = 0.01f;           // small value 

// smoothing kernels defined in Müller and their gradients
// adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
const static float POLY6 = 4.f / (M_PI * pow(H, 8.f));
const static float SPIKY_GRAD = -10.f / (M_PI * pow(H, 5.f));
const static float VISC_LAP = 40.f / (M_PI * pow(H, 5.f));

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;

// particle data structure
// stores position, velocity, and force for integration, New : also define the normal for each particle
// stores density (rho) and pressure values for SPH
struct Particle
{
	Particle(float _x, float _y, float _z) : x(_x, _y, _z), v(0.f, 0.f, 0.f), f(0.f, 0.f, 0.f), n(0.f, 0.f, 0.f), rho(0), p(0.f) {}
	Vector3d x, v, f, n;
	float rho, p;
	Vector3d color; // for rendering
	bool onSurface = false;
};

// solver data
static vector<Particle> particles;
const static int DAM_PARTICLES = 2500;

// rendering projection parameters
const static int WINDOW_WIDTH = 2000;
const static int WINDOW_HEIGHT = 1000;
const static int WINDOW_DEPTH = 900;
const static double VIEW_WIDTH = 1.5 * 800.f;
const static double VIEW_HEIGHT = 1.5 * 600.f;
const static double VIEW_DEPTH = 1.5 * 600.f;

//new create thing
GLFWwindow* window;


void InitSPH(void)
{
	cout << "initializing dam break with " << DAM_PARTICLES << " particles" << endl;
	for (float z = EPS; z < VIEW_DEPTH - EPS * 2.f; z += H)
	{
		for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H)
		{
			for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H)
			{
				if (particles.size() < DAM_PARTICLES)
				{
					float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
					particles.push_back(Particle(x + jitter, y+jitter, z));
				}
				else
				{
					return;
				}
			}
		}
	}
}

void Integrate(void)
{
	for (auto& p : particles)
	{
		// forward Euler integration
		p.v += DT * p.f / p.rho;
		p.x += DT * p.v;

		// enforce boundary conditions
		if (p.x(0) - EPS < 0.f)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = EPS;
		}
		if (p.x(0) + EPS > VIEW_WIDTH)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = VIEW_WIDTH - EPS;
		}
		if (p.x(1) - EPS < 0.f)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = EPS;
		}
		if (p.x(1) + EPS > VIEW_HEIGHT)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = VIEW_HEIGHT - EPS;
		}
		if (p.x(2) - EPS < 0.f)
		{
			p.v(2) *= BOUND_DAMPING;
			p.x(2) = EPS;
		}
		if (p.x(2) + EPS > VIEW_DEPTH)
		{
			p.v(2) *= BOUND_DAMPING;
			p.x(2) = VIEW_DEPTH - EPS;
		}
	}
}


void ComputeDensityPressure(void)
{
	for (auto& pi : particles)
	{
		pi.rho = 0.f;
		for (auto& pj : particles)
		{
			Vector3d rij = pj.x - pi.x;
			float r2 = rij.squaredNorm();

			if (r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
			}
		}
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}


float kernel(float r, float h) {
	if (r > h) {
		return 0.0f;
	}
	float q = h-r;
	float s3 = 1.0f / ( M_PI * h * h * h);
	if (q < 1.0f) {
		float q2 = q * q;
		float q3 = q2 * q;
		return s3 * (1.0f - 1.5f * q2 + 0.5f * q3);
	}
	else if (q < 2.0f) {
		float q2 = q * q;
		float q3 = q2 * q;
		return s3 * 0.25f * powf(2.0f - q, 3.0f);
	}
	else {
		return 0.0f;
	}
}

float gradient_kernel(float r) {
	float r1 = r + epsilon;
	float r2 = r - epsilon;
	float g_k = (kernel(r1, H) - kernel(r2, H)) / (2 * epsilon);
	return g_k;
}

void ComputeNormal(void) {
	Vector3d normal(0.f, 0.f, 0.f);
	float r;
	Vector3d rij(0.f, 0.f, 0.f);
	for (auto& pi : particles)
	{
		for (auto& pj : particles)
		{
			rij = pj.x - pi.x;
			r = rij.norm();
			if (r < H) {
				float g = gradient_kernel(r);
				normal += MASS * (1.f / pj.rho)* g * pj.n;
			}
		}
		if (normal.norm() > 0.3f) { //threshold l for optimazation. i.e. only consider the particle on the surface
			pi.onSurface = true;
		}

		pi.n = H * normal;
	}
	
}

void ComputeForces(void)
{
	for (auto& pi : particles)
	{
		Vector3d fpress(0.f, 0.f, 0.f);
		Vector3d fvisc(0.f, 0.f, 0.f);
		Vector3d fsurf(0.f, 0.f, 0.f);
		for (auto& pj : particles)
		{
			if (&pi == &pj) continue;
			Vector3d rij = pj.x - pi.x;
			float r = rij.norm();
			if (r < H)
			{
				// compute pressure force contribution
				fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 3.f);
				// compute viscosity force contribution
				fvisc += VISC * MASS * (pj.v - pi.v) / pj.rho * VISC_LAP * (H - r);
				// compute the surface tension
				if (pi.onSurface) fsurf -= sigma * MASS *(pi.n-pj.n);
			}
		}

		// add the external force, here is the gravity
		Vector3d fgrav = G * MASS / pi.rho;
		pi.f = fpress + fvisc + fgrav + fsurf;
	}
}

void ComputeColor() {
	// normalize all the velocity between 0 and 1
	float max_v = 0.f;
	for (auto& pi : particles) {
		//float v = std::sqrt(pi.x[0] * pi.x[0] + pi.x[1] * pi.x[1]); // position space as the texture space
		float v = std::sqrt(pi.v[0] * pi.v[0] + pi.v[1] * pi.v[1]); // velocity space as the texture space
		max_v = std::max(max_v, v);
	}

	//convert velocity(or position) to HSV, then convert HSV to RGB
	for (auto& pi : particles) {
		//float v = std::sqrt(pi.x[0] * pi.x[0] + pi.x[1] * pi.x[1]);
		float v = std::sqrt(pi.v[0] * pi.v[0] + pi.v[1] * pi.v[1]);
		
		if (v > 0) {
			// convert velocity to hue
			// create a linear mapping  : (x,y,z) -> (x,y), then another mapping : Cartesian coordinate -> polar coordinate to reperesnt hue as angle
			// angle in range [0, 2*pi]
			//float theta = std::atan2(pi.x[1], pi.x[0]) + M_PI; // position space as the texture space
			float theta = std::atan2(pi.v[1], pi.v[0]) + M_PI; // velocity space as the texture space
			float h = theta / (2 * M_PI); // then convert into [0,1] 
			// convert velocity to satuaration, in range [0,1]
			float s = v/max_v;
			// set the value to be a constant. here we use full brightness (value = 1), for prettier animation
			float value = 1.f;

			// Convert HSV to RGB
			h = h * 6.f;
			float m = value - s;
			float f = std::abs(std::fmod(h, 2.f) - 1.f); //make sure f < 1 so that we won't get negative value for RGB
			float k = value - s * (1.f - f);
			float n = value - s * f;
			if (h < 1.0) {
				pi.color = Vector3d(value, k, m);
			}
			else if (h < 2.0) {
				pi.color = Vector3d(n, value, m);
			}
			else if (h < 3.0) {
				pi.color = Vector3d(m, value, k);
			}
			else if (h < 4.0) {
				pi.color = Vector3d(m, n, value);
			}
			else if (h < 5.0) {
				pi.color = Vector3d(k, m, value);
			}
			else {
				pi.color = Vector3d(value, m, n);
			}
		}
		else {
			// if v = 0, then the color is black
			pi.color = Vector3d(0.0f, 0.0f, 0.0f);
		}
	}
}


void Update(void)
{
	ComputeDensityPressure();
	ComputeNormal();
	ComputeForces();
	Integrate();
	ComputeColor();

	glfwPostEmptyEvent();
}


void InitGL(void)
{
	glClearColor(0.9f, 0.9f, 0.9f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(H/2.f);

	//viewport
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	float ratio = width / (float)(height);
	//projection
	glm::mat4 projection_m = glm::perspective(80.f, ratio, 0.1f, 1000.0f);
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projection_m));

	//modeling and viewing
	glm::mat4 modelview_m = glm::lookAt(glm::vec3(0, 0, 700), glm::vec3(400, 0, 0), glm::vec3(0,-20, 0));
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(modelview_m));
}

void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBegin(GL_POINTS);

	for (auto& p : particles)
	{
		glColor4f(p.color[0], p.color[1], p.color[2], 1);
		glVertex3f(p.x(0), p.x(1), p.x(2));
	}
	glEnd();

	// swap buffers
	glfwSwapBuffers(window);
}

void Keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action != GLFW_PRESS) return;
	if (key == GLFW_KEY_R){
		particles.clear();
		InitSPH();
	}
}

int main(int argc, char** argv)
{
	if (!glfwInit()) return -1;

	// create a window
	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "candy ocean", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);

	// keyboard function
	glfwSetKeyCallback(window, Keyboard);

	// initialize OpenGL 
	InitGL();

	// initialize SPH system
	InitSPH();

	// loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		// render
		Render(); // swap buffer is set inside the render()

		//update particles
		Update();

		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}


