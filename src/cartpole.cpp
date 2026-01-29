#include "cartpole.hpp"
#include <cmath>

namespace cartpole {

State derivatives(const State& s, double u, const Params& p) {
    const double M = p.m_cart;
    const double m = p.m_pole;
    const double L = p.L;
    const double g = p.g;
    const double mu_c = p.mu_cart;
    const double mu_p = p.mu_pole;
    
    const double theta = s.theta;
    const double theta_dot = s.theta_dot;
    const double x_dot = s.x_dot;
    
    const double sin_t = std::sin(theta);
    const double cos_t = std::cos(theta);
    const double cos_t_sq = cos_t * cos_t;
    
    const double denom = M + m - m * cos_t_sq; 

    const double friction_cart = mu_c * std::tanh(1000.0 * x_dot);
    
    const double friction_pole = mu_p * theta_dot;

    double theta_ddot = (
        (M + m) * g * sin_t
        - cos_t * (u + m * L * theta_dot * theta_dot * sin_t - friction_cart)
        - friction_pole * (M + m) / (m * L)
    ) / (L * denom);
    

    double x_ddot = (
        u 
        + m * L * theta_dot * theta_dot * sin_t
        - m * L * theta_ddot * cos_t
        - friction_cart
    ) / (M + m);
    

    x_ddot = (
        u 
        + m * L * theta_dot * theta_dot * sin_t
        - m * g * sin_t * cos_t
        - friction_cart
        + friction_pole * cos_t / L
    ) / denom;
    
    theta_ddot = (
        (M + m) * g * sin_t
        - m * L * theta_dot * theta_dot * sin_t * cos_t
        - (u - friction_cart) * cos_t
        - friction_pole * (M + m) / (m * L)
    ) / (L * denom);
    
    State deriv;
    deriv.x = x_dot;            
    deriv.x_dot = x_ddot;     
    deriv.theta = theta_dot;   
    deriv.theta_dot = theta_ddot; 
    
    return deriv;
}

double normalize_angle(double angle) {
    while (angle > M_PI) {
        angle -= 2.0 * M_PI;
    }
    while (angle < -M_PI) {
        angle += 2.0 * M_PI;
    }
    return angle;
}

}