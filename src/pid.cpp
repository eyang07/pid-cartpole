#include "pid.hpp"
#include <algorithm>
#include <cmath>

namespace pid {

PIDController::PIDController()
    : config_()
    , integral_(0.0)
    , prev_measurement_(0.0)
    , first_update_(true)
    , last_error_(0.0)
    , derivative_(0.0)
    , p_term_(0.0)
    , i_term_(0.0)
    , d_term_(0.0)
{}

PIDController::PIDController(const PIDConfig& config)
    : config_(config)
    , integral_(0.0)
    , prev_measurement_(0.0)
    , first_update_(true)
    , last_error_(0.0)
    , derivative_(0.0)
    , p_term_(0.0)
    , i_term_(0.0)
    , d_term_(0.0)
{}

double PIDController::update(double measurement, double dt) {
    double error = config_.setpoint - measurement;
    last_error_ = error;
    
    p_term_ = config_.Kp * error;
    
    integral_ += error * dt;
    integral_ = std::clamp(integral_, config_.integral_min, config_.integral_max);
    i_term_ = config_.Ki * integral_;
    
    if (first_update_) {
        derivative_ = 0.0;
        first_update_ = false;
    } else {
        derivative_ = -(measurement - prev_measurement_) / dt;
    }
    prev_measurement_ = measurement;
    d_term_ = config_.Kd * derivative_;
    
    double output = p_term_ + i_term_ + d_term_;
    
    output = std::clamp(output, config_.output_min, config_.output_max);
    
    return output;
}

void PIDController::reset() {
    integral_ = 0.0;
    prev_measurement_ = 0.0;
    first_update_ = true;
    last_error_ = 0.0;
    derivative_ = 0.0;
    p_term_ = 0.0;
    i_term_ = 0.0;
    d_term_ = 0.0;
}

void PIDController::set_config(const PIDConfig& config) {
    config_ = config;
}

CascadedPID::CascadedPID()
    : angle_pid_()
    , position_pid_()
{}

CascadedPID::CascadedPID(const PIDConfig& angle_config, const PIDConfig& position_config)
    : angle_pid_(angle_config)
    , position_pid_(position_config)
{}

double CascadedPID::update(double theta, double x, double dt) {
    double angle_offset = position_pid_.update(x, dt);
    
    const double max_angle_offset = 0.15; 
    angle_offset = std::clamp(angle_offset, -max_angle_offset, max_angle_offset);
    
    double target_angle = angle_pid_.get_setpoint() + angle_offset;
    
    double original_setpoint = angle_pid_.get_setpoint();
    angle_pid_.set_setpoint(target_angle);
    double force = angle_pid_.update(theta, dt);
    angle_pid_.set_setpoint(original_setpoint);
    
    return force;
}

void CascadedPID::reset() {
    angle_pid_.reset();
    position_pid_.reset();
}

} // namespace pid