//
//  TimeGenerator.h
//  CTRW
//
//  Created by Tomas Aquino on 10/23/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

//	Generic waiting time generators

#ifndef TimeGenerator_h
#define TimeGenerator_h

namespace ctrw
{
  template <typename val_type>
  class TimeGenerator_Step
  {
    //  Deterministic time step
    val_type dt;

  public:
    using value_type = val_type;

    TimeGenerator_Step(val_type dt = 0.)
    : dt(dt)
    {}

    void time_step(val_type dt)
    { this->dt = dt; }

    val_type time_step() const
    { return dt; }

    template <typename TT>
    val_type operator() (TT const& = {})
    { return dt; }
  };
  
  template <typename Reference_angle>
  class TimeGenerator_AngleClasses_Exponential_1d
  {
  public:
    const std::vector<double> angle_classes;
    const std::vector<double> mean_times;
    
    TimeGenerator_AngleClasses_Exponential_1d
    (std::vector<double> angle_classes, std::vector<double> mean_times,
     Reference_angle const& reference_angle)
    : angle_classes{ angle_classes }
    , mean_times{ mean_times }
    , reference_angle{ reference_angle }
    {}
    
    template <typename State>
    double operator()(State const& state)
    {
      std::size_t angle_idx = angle_class(state.angle-reference_angle(state));
      
      return mean_times[angle_idx]*dist(rng);
    }
    
    std::size_t angle_class(double angle)
    {
      return std::lower_bound(angle_classes.begin(), angle_classes.end(), angle) -
        angle_classes.begin();
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution<double> dist{ 1. };
    Reference_angle const& reference_angle;
  };
}



#endif /* TimeGenerator_h */
