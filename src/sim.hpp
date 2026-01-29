#ifndef SIM_HPP
#define SIM_HPP

#include "cartpole.hpp"
#include <functional>

namespace sim {

using ControlFn = std::function<double(const cartpole::State&, double)>;

cartpole::State rk4_step(
    const cartpole::State& s,
    double u,
    const cartpole::Params& p,
    double dt
);

struct SimConfig {
    double dt;           
    double t_end;    
    double u_max;  
    bool apply_disturbance; 
    double disturbance_time;
    double disturbance_mag; 
    
    SimConfig()
        : dt(0.002)
        , t_end(10.0)
        , u_max(100.0)
        , apply_disturbance(false)
        , disturbance_time(2.0)
        , disturbance_mag(5.0)
    {}
};

struct SimStep {
    double t;             
    cartpole::State state;
    double u;       
    double u_raw;
};

std::vector<SimStep> run_simulation(
    const cartpole::State& initial,
    const cartpole::Params& params,
    const SimConfig& config,
    ControlFn control
);

double saturate(double value, double limit);

}

#endif // SIM_HPP