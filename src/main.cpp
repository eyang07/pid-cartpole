#include "cartpole.hpp"
#include "sim.hpp"
#include "pid.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

void write_csv(const std::string& filename, const std::vector<sim::SimStep>& results) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing\n";
        return;
    }
    
    file << "time,x,x_dot,theta,theta_dot,theta_deg,u,u_raw\n";
    
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < results.size(); i += 5) {
        const auto& r = results[i];
        file << r.t << ","
             << r.state.x << ","
             << r.state.x_dot << ","
             << r.state.theta << ","
             << r.state.theta_dot << ","
             << r.state.theta * 180.0 / M_PI << ","
             << r.u << ","
             << r.u_raw << "\n";
    }
    
    file.close();
}

struct Metrics {
    bool success;         
    double max_angle_deg;
    double settling_time;
    double mean_control;
    double max_control;
    double final_angle_deg;
    double final_x;
};

Metrics compute_metrics(const std::vector<sim::SimStep>& results, double sim_duration) {
    Metrics m;
    m.success = true;
    m.max_angle_deg = 0.0;
    m.settling_time = -1.0;
    m.mean_control = 0.0;
    m.max_control = 0.0;
    
    const double success_threshold_deg = 12.0;
    const double settling_threshold_deg = 2.0;
    
    double control_sum = 0.0;
    bool settled = false;
    double settle_start = -1.0;
    
    for (const auto& r : results) {
        double angle_deg = std::abs(r.state.theta) * 180.0 / M_PI;
        
        if (angle_deg > success_threshold_deg) {
            m.success = false;
        }
        
        m.max_angle_deg = std::max(m.max_angle_deg, angle_deg);
        
        if (angle_deg <= settling_threshold_deg) {
            if (settle_start < 0) {
                settle_start = r.t;
            }
        } else {
            settle_start = -1.0;
        }
        
        double abs_u = std::abs(r.u);
        control_sum += abs_u;
        m.max_control = std::max(m.max_control, abs_u);
    }
    
    if (settle_start >= 0) {
        m.settling_time = settle_start;
    }
    
    m.mean_control = control_sum / results.size();
    m.final_angle_deg = results.back().state.theta * 180.0 / M_PI;
    m.final_x = results.back().state.x;
    
    return m;
}

void print_metrics(const Metrics& m, const std::string& label = "") {
    if (!label.empty()) {
        std::cout << label << ":\n";
    }
    std::cout << "  Success (within ±12°): " << (m.success ? "YES" : "NO") << "\n"
              << "  Max angle:    " << std::fixed << std::setprecision(2) << m.max_angle_deg << "°\n"
              << "  Settling time: ";
    if (m.settling_time >= 0) {
        std::cout << m.settling_time << " s\n";
    } else {
        std::cout << "never settled\n";
    }
    std::cout << "  Final angle:  " << m.final_angle_deg << "°\n"
              << "  Final x:      " << m.final_x << " m\n"
              << "  Mean |u|:     " << m.mean_control << " N\n"
              << "  Max |u|:      " << m.max_control << " N\n";
}

int main(int argc, char* argv[]) {
    using namespace cartpole;
    using namespace sim;
    using namespace pid;
    
    fs::create_directories("data");
    
    Params params;
    std::cout << "System Parameters:\n"
              << "  Cart mass:    " << params.m_cart << " kg\n"
              << "  Pole mass:    " << params.m_pole << " kg\n"
              << "  Pole length:  " << params.L * 2 << " m (half-length: " << params.L << " m)\n"
              << "  Gravity:      " << params.g << " m/s²\n\n";
    
    SimConfig sim_cfg;
    sim_cfg.dt = 0.002;    
    sim_cfg.t_end = 10.0;  
    sim_cfg.u_max = 100.0; 
    
    std::cout << "Simulation Config:\n"
              << "  Timestep:     " << sim_cfg.dt * 1000 << " ms\n"
              << "  Duration:     " << sim_cfg.t_end << " s\n"
              << "  Force limit:  ±" << sim_cfg.u_max << " N\n\n";
    
    PIDConfig angle_cfg(-200.0, -2.0, -30.0);
    angle_cfg.setpoint = 0.0;
    angle_cfg.output_min = -sim_cfg.u_max;
    angle_cfg.output_max = sim_cfg.u_max;
    angle_cfg.integral_min = -20.0;
    angle_cfg.integral_max = 20.0;
    
    PIDConfig pos_cfg(0.02, 0.002, 0.08);
    pos_cfg.setpoint = 0.0;
    pos_cfg.output_min = -0.2;
    pos_cfg.output_max = 0.2;
    pos_cfg.integral_min = -0.5;
    pos_cfg.integral_max = 0.5;
    
    std::cout << "PID Gains:\n"
              << "  Angle:    Kp=" << angle_cfg.Kp << ", Ki=" << angle_cfg.Ki 
              << ", Kd=" << angle_cfg.Kd << "\n"
              << "  Position: Kp=" << pos_cfg.Kp << ", Ki=" << pos_cfg.Ki 
              << ", Kd=" << pos_cfg.Kd << "\n\n";
    
    {
        CascadedPID controller(angle_cfg, pos_cfg);
        State initial(0, 0, 10.0 * M_PI / 180.0, 0); 
        
        auto results = run_simulation(initial, params, sim_cfg,
            [&controller, &sim_cfg](const State& s, double) {
                return controller.update(s.theta, s.x, sim_cfg.dt);
            }
        );
        
        write_csv("data/demo_run.csv", results);
        
        Metrics m = compute_metrics(results, sim_cfg.t_end);
        print_metrics(m);
        std::cout << "  Saved to: data/demo_run.csv\n\n";
    }
    
    {
        CascadedPID controller(angle_cfg, pos_cfg);
        State initial(0, 0, 5.0 * M_PI / 180.0, 0);
        
        SimConfig dist_cfg = sim_cfg;
        dist_cfg.apply_disturbance = true;
        dist_cfg.disturbance_time = 2.0;
        dist_cfg.disturbance_mag = 3.0;
        
        auto results = run_simulation(initial, params, dist_cfg,
            [&controller, &sim_cfg](const State& s, double) {
                return controller.update(s.theta, s.x, sim_cfg.dt);
            }
        );
        
        write_csv("data/demo_disturbance.csv", results);
        
        Metrics m = compute_metrics(results, sim_cfg.t_end);
        print_metrics(m);
        std::cout << "  Saved to: data/demo_disturbance.csv\n\n";
    }
    
    std::cout << "Randomized Test Suite (100 trials)\n";
    std::cout << "  Initial angle:    [-15°, 15°]\n";
    std::cout << "  Initial ang. vel: [-0.5, 0.5] rad/s\n";
    std::cout << "  50% with disturbance at t=2s\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> angle_dist(-15.0, 15.0);
    std::uniform_real_distribution<> angvel_dist(-0.5, 0.5);
    std::bernoulli_distribution disturbance_dist(0.5);
    
    const int num_trials = 100;
    int successes = 0;
    std::vector<double> settling_times;
    std::vector<double> max_angles;
    std::vector<double> mean_controls;
    
    std::ofstream metrics_file("data/trial_metrics.csv");
    metrics_file << "trial,theta0_deg,thetadot0,disturbance,success,max_angle_deg,"
                 << "settling_time,final_angle_deg,final_x,mean_u,max_u\n";
    
    for (int trial = 0; trial < num_trials; trial++) {
        double theta0_deg = angle_dist(gen);
        double theta0 = theta0_deg * M_PI / 180.0;
        double thetadot0 = angvel_dist(gen);
        bool apply_dist = disturbance_dist(gen);
        
        CascadedPID controller(angle_cfg, pos_cfg);
        State initial(0, 0, theta0, thetadot0);
        
        SimConfig trial_cfg = sim_cfg;
        trial_cfg.apply_disturbance = apply_dist;
        trial_cfg.disturbance_time = 2.0;
        trial_cfg.disturbance_mag = 3.0;
        
        auto results = run_simulation(initial, params, trial_cfg,
            [&controller, &sim_cfg](const State& s, double) {
                return controller.update(s.theta, s.x, sim_cfg.dt);
            }
        );
        
        Metrics m = compute_metrics(results, sim_cfg.t_end);
        
        if (m.success) successes++;
        if (m.settling_time >= 0) settling_times.push_back(m.settling_time);
        max_angles.push_back(m.max_angle_deg);
        mean_controls.push_back(m.mean_control);
        
        metrics_file << trial << ","
                     << theta0_deg << ","
                     << thetadot0 << ","
                     << (apply_dist ? 1 : 0) << ","
                     << (m.success ? 1 : 0) << ","
                     << m.max_angle_deg << ","
                     << m.settling_time << ","
                     << m.final_angle_deg << ","
                     << m.final_x << ","
                     << m.mean_control << ","
                     << m.max_control << "\n";
        
        if (trial < 5) {
            write_csv("data/trial_" + std::to_string(trial) + ".csv", results);
        }
    }
    
    metrics_file.close();
    
    double success_rate = 100.0 * successes / num_trials;
    
    double mean_settling = 0;
    if (!settling_times.empty()) {
        for (double t : settling_times) mean_settling += t;
        mean_settling /= settling_times.size();
    }
    
    double mean_max_angle = 0;
    for (double a : max_angles) mean_max_angle += a;
    mean_max_angle /= max_angles.size();
    
    double mean_mean_control = 0;
    for (double u : mean_controls) mean_mean_control += u;
    mean_mean_control /= mean_controls.size();
    
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "\nResults:\n"
              << "  Success rate:      " << success_rate << "% (" << successes << "/" << num_trials << ")\n"
              << "  Mean settling time: " << mean_settling << " s (of " << settling_times.size() << " that settled)\n"
              << "  Mean max angle:    " << mean_max_angle << "°\n"
              << "  Mean control effort: " << mean_mean_control << " N\n"
              << "\n  Saved to: data/trial_metrics.csv\n";
    ;
    std::cout << "Edge Case Tests\n";
    
    {
        CascadedPID controller(angle_cfg, pos_cfg);
        State initial(0, 0, 14.0 * M_PI / 180.0, 0);
        
        auto results = run_simulation(initial, params, sim_cfg,
            [&controller, &sim_cfg](const State& s, double) {
                return controller.update(s.theta, s.x, sim_cfg.dt);
            }
        );
        
        write_csv("data/edge_large_angle.csv", results);
        Metrics m = compute_metrics(results, sim_cfg.t_end);
        print_metrics(m, "Large initial angle (14°)");
    }
    
    {
        CascadedPID controller(angle_cfg, pos_cfg);
        State initial(0, 0, 0.0, 1.0);
        
        auto results = run_simulation(initial, params, sim_cfg,
            [&controller, &sim_cfg](const State& s, double) {
                return controller.update(s.theta, s.x, sim_cfg.dt);
            }
        );
        
        write_csv("data/edge_angular_velocity.csv", results);
        Metrics m = compute_metrics(results, sim_cfg.t_end);
        print_metrics(m, "\nInitial angular velocity (1 rad/s)");
    }
    
    return 0;
}