#ifndef PID_HPP
#define PID_HPP

namespace pid {

struct PIDConfig {
    double Kp;
    double Ki;          
    double Kd;
    double setpoint;
    double output_min; 
    double output_max;   
    double integral_min;  
    double integral_max;  
    
    PIDConfig()
        : Kp(0.0)
        , Ki(0.0)
        , Kd(0.0)
        , setpoint(0.0)
        , output_min(-100.0)
        , output_max(100.0)
        , integral_min(-50.0)
        , integral_max(50.0)
    {}
    
    PIDConfig(double kp, double ki, double kd,
              double sp = 0.0,
              double out_min = -100.0, double out_max = 100.0,
              double int_min = -50.0, double int_max = 50.0)
        : Kp(kp)
        , Ki(ki)
        , Kd(kd)
        , setpoint(sp)
        , output_min(out_min)
        , output_max(out_max)
        , integral_min(int_min)
        , integral_max(int_max)
    {}
};

class PIDController {
public:
    PIDController();
    explicit PIDController(const PIDConfig& config);
    

    double update(double measurement, double dt);
    
    void reset();
    
    void set_config(const PIDConfig& config);
    const PIDConfig& get_config() const { return config_; }
    
    void set_setpoint(double sp) { config_.setpoint = sp; }
    double get_setpoint() const { return config_.setpoint; }
    
    double get_error() const { return last_error_; }
    double get_integral() const { return integral_; }
    double get_derivative() const { return derivative_; }
    double get_p_term() const { return p_term_; }
    double get_i_term() const { return i_term_; }
    double get_d_term() const { return d_term_; }
    
private:
    PIDConfig config_;
    
    double integral_;          
    double prev_measurement_;  
    bool first_update_;    
    
    double last_error_;
    double derivative_;
    double p_term_;
    double i_term_;
    double d_term_;
};

class CascadedPID {
public:
    CascadedPID();
    CascadedPID(const PIDConfig& angle_config, const PIDConfig& position_config);
    
    double update(double theta, double x, double dt);
    
    void reset();
    
    PIDController& angle_controller() { return angle_pid_; }
    PIDController& position_controller() { return position_pid_; }
    
private:
    PIDController angle_pid_;    
    PIDController position_pid_;   
};

}

#endif // PID_HPP