#ifndef CARTPOLE_HPP
#define CARTPOLE_HPP

namespace cartpole {

// State vector: [x, x_dot, theta, theta_dot]
// x: cart position (m)
// x_dot: cart velocity (m/s)
// theta: pole angle from vertical (rad), 0 = upright
// theta_dot: pole angular velocity (rad/s)
struct State {
    double x;
    double x_dot;
    double theta;
    double theta_dot;
    
    State() : x(0), x_dot(0), theta(0), theta_dot(0) {}
    State(double x_, double x_dot_, double theta_, double theta_dot_)
        : x(x_), x_dot(x_dot_), theta(theta_), theta_dot(theta_dot_) {}
};

// System parameters
struct Params {
    double m_cart;      // cart mass (kg)
    double m_pole;      // pole mass (kg)
    double L;           // pole half-length (m) - distance to center of mass
    double g;           // gravitational acceleration (m/s^2)
    double mu_cart;     // cart friction coefficient
    double mu_pole;     // pole friction coefficient (at pivot)
    
    Params()
        : m_cart(1.0)
        , m_pole(0.1)
        , L(0.5)
        , g(9.81)
        , mu_cart(0.0005)
        , mu_pole(0.000002)
    {}
    
    Params(double mc, double mp, double l, double gravity,
           double mu_c = 0.0005, double mu_p = 0.000002)
        : m_cart(mc)
        , m_pole(mp)
        , L(l)
        , g(gravity)
        , mu_cart(mu_c)
        , mu_pole(mu_p)
    {}
    
    double total_mass() const { return m_cart + m_pole; }
};

State derivatives(const State& s, double u, const Params& p);

double normalize_angle(double angle);

} 

#endif // CARTPOLE_HPP