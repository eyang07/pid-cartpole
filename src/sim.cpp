#include "sim.hpp"
#include <cmath>
#include <algorithm>

namespace sim {

static cartpole::State state_add(const cartpole::State& a, const cartpole::State& b) {
    return cartpole::State(
        a.x + b.x,
        a.x_dot + b.x_dot,
        a.theta + b.theta,
        a.theta_dot + b.theta_dot
    );
}

static cartpole::State state_scale(const cartpole::State& s, double k) {
    return cartpole::State(
        s.x * k,
        s.x_dot * k,
        s.theta * k,
        s.theta_dot * k
    );
}

// RK4 integration step
// 
//   k1 = f(t, y)
//   k2 = f(t + dt/2, y + dt/2 * k1)
//   k3 = f(t + dt/2, y + dt/2 * k2)
//   k4 = f(t + dt, y + dt * k3)
//   y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
//
// For our system, f computes the state derivatives given state and control.
// We hold the control force constant over the timestep (zero-order hold).

cartpole::State rk4_step(
    const cartpole::State& s,
    double u,
    const cartpole::Params& p,
    double dt
) {
    cartpole::State k1 = cartpole::derivatives(s, u, p);
    
    cartpole::State s2 = state_add(s, state_scale(k1, dt / 2.0));
    cartpole::State k2 = cartpole::derivatives(s2, u, p);
    
    cartpole::State s3 = state_add(s, state_scale(k2, dt / 2.0));
    cartpole::State k3 = cartpole::derivatives(s3, u, p);
    
    cartpole::State s4 = state_add(s, state_scale(k3, dt));
    cartpole::State k4 = cartpole::derivatives(s4, u, p);
    
    cartpole::State weighted_sum(
        k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x,
        k1.x_dot + 2.0 * k2.x_dot + 2.0 * k3.x_dot + k4.x_dot,
        k1.theta + 2.0 * k2.theta + 2.0 * k3.theta + k4.theta,
        k1.theta_dot + 2.0 * k2.theta_dot + 2.0 * k3.theta_dot + k4.theta_dot
    );
    
    cartpole::State new_state = state_add(s, state_scale(weighted_sum, dt / 6.0));
    
    new_state.theta = cartpole::normalize_angle(new_state.theta);
    
    return new_state;
}

double saturate(double value, double limit) {
    return std::clamp(value, -limit, limit);
}

std::vector<SimStep> run_simulation(
    const cartpole::State& initial,
    const cartpole::Params& params,
    const SimConfig& config,
    ControlFn control
) {
    size_t n_steps = static_cast<size_t>(config.t_end / config.dt) + 1;
    std::vector<SimStep> results;
    results.reserve(n_steps);
    
    cartpole::State state = initial;
    double t = 0.0;
    
    while (t <= config.t_end) {
        double u_raw = control(state, t);
        
        double u = saturate(u_raw, config.u_max);
        
        if (config.apply_disturbance) {
            double t_dist = config.disturbance_time;
            if (t >= t_dist && t < t_dist + config.dt) {
                u += config.disturbance_mag / config.dt;
            }
        }
        
        SimStep step;
        step.t = t;
        step.state = state;
        step.u = u;
        step.u_raw = u_raw;
        results.push_back(step);
        
        state = rk4_step(state, u, params, config.dt);
        t += config.dt;
    }
    
    return results;
}

}