//
//  Measurer.h
//  CTRW
//
//  Created by Tomas Aquino on 14/1/20.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

#ifndef Measurer_CTRW_h
#define Measurer_CTRW_h

#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <unordered_set>
#include "general/Operations.h"
#include "general/useful.h"
#include "Stochastic/CTRW/StateGetter.h"

namespace ctrw
{
  // Measurer_Name objects implement operator() to
  // measure and output at the same time
  // Measurer_Store_Name objects implement update method to
  // measure and store, get method to get stored values,
  // and print method to output stored values
  
  // Output based on particle states
  // Output objects in methods must take
  // a state, an std::ofstream object, and a delimiter
  // to output wanted state information
  class Measurer_State
  {
  public:
    struct New{};   // To pick measure new state
    struct Old{};   // To pick measuring old state
    struct Both{};  // To pick measuring both states
    
    // Construct given output filename,
    // output scientific notation precision,
    // and output file delimiter
    Measurer_State
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
    }
    
    // Destructor cleans up output object
    ~Measurer_State()
    { output.close(); }
      
    // Output particle new states
    template <typename Subject, typename OutputState>
    void operator()
    (Subject const& subject, OutputState output_state, New)
    {
      bool delim = 0;
      for (auto const& part : subject.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_new(), output, delimiter);
        delim = 1;
      }
      output << "\n";
    }
      
    // Output particle old states
    template <typename Subject, typename OutputState>
    void operator()
    (Subject const& subject, OutputState output_state, Old)
    {
      bool delim = 0;
      for (auto const& part : subject.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_old(), output, delimiter);
        delim = 1;
      }
      output << "\n";
    }
    
    // Output particle new and old states
    template <typename Subject, typename OutputState>
    void operator()
    (Subject const& subject, OutputState output_state, Both)
    {
      bool delim = 0;
      for (auto const& part : subject.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_new(), output, delimiter);
        output << delimiter;
        output_state(part.state_old(), output, delimiter);
      }
      output << "\n";
    }
      
    // Output each particle state information preceded by scalar value
    template <typename Subject, typename OutputState,
    typename Value,typename Which>
    void operator()
    (Subject const& subject, OutputState output_state, Value tag,
     Which which)
    {
      output << tag << delimiter;
      (*this)(subject, output_state, which);
    }
    
  private:
    std::ofstream output;   // Handle output to file
    std::string delimiter;  // Delimiter between state quantities in each file line
  };
      
  // Output sum of wanted quantity over all particles at each measurement
  // Getter objects in methods must take a particle and output wanted quantity
  // Wanted quantity must be summable by useful::plus_InPlace
  // and printable by useful::print
  // Subject must implement cbegin() and size()
  class Measurer_Total
  {
  public:
    // Construct given output filename,
    // output scientific notation precision,
    // and output file delimiter
    Measurer_Total
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
    }
    
    // Destructor cleans up output object
    ~Measurer_Total()
    { output.close(); }
    
    // Output total
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
        return;
      
      auto val = get(*subject.cbegin());
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend(); ++part_it)
        operation::plus_InPlace(val, get(*part_it));
      useful::print(output, val, 0, delimiter);
      output << "\n";
    }
    
    // Output total preceded by scalar value
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
    private:
      std::ofstream output;   // Handle output to file
      std::string delimiter;  // Delimiter between state quantities in each file line
  };
      
  // Output arithmetic mean of wanted quantity over all particles at each measurement
  // Getter objects in methods must take a particle and output wanted quantity
  // Wanted quantity must be summable by useful::plus_InPlace,
  // divisible by scalar by useful::div_scalar_InPlace,
  // and printable by useful::print
  // Subject must implement cbegin() and size()
  class Measurer_Mean
  {
  public:
    // Construct given output filename,
    // output scientific notation precision,
    // and output file delimiter
    Measurer_Mean
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
    }
    
    // Destructor cleans up output object
    ~Measurer_Mean()
    { output.close(); }
    
    // Output mean
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
        return;
      
      auto val = get(*subject.cbegin());
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend(); ++part_it)
        operation::plus_InPlace(val, get(*part_it));
      operation::div_scalar_InPlace(val, subject.size());
      useful::print(output, val, 0, delimiter);
      output << "\n";
    }
    
    // Output mean preceded by scalar value
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
    private:
      std::ofstream output;   // Handle output to file
      std::string delimiter;  // Delimiter between state quantities in each file line
  };
  
  // Output based on particles
  // Getter objects in methods must take a particle
  // and return wanted quantity
  // Wanted quantity must be printable by useful::print
  // Subject must implement particles()
  class Measurer_Particle
  {
  public:
    // Construct given output filename,
    // output scientific notation precision,
    // and output file delimiter
    Measurer_Particle
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
    }
    
    // Destructor cleans up output object
    ~Measurer_Particle()
    { output.close(); }
    
    // Output each particle's info
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      bool delim = 0;
      for (auto const& part : subject.particles())
      {
        useful::print(output, get(part), delim, delimiter);
        delim = 1;
      }
      output << "\n";
    }
    
    // Output each particle's info preceded by scalar value
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
  private:
    std::ofstream output;   // Handle output to file
    std::string delimiter;  // Delimiter between state quantities in each file line
  };
        
  // For each element ii in subject collection,
  // Iterates starting from subject.cbegin(ii) to subject.cend(ii)
  // and outputs size(ii) followed by corresponding values
  // Subject must implement size(), size(std::size_t), cbegin(std::size_t), and cend(std::size_t)
  // Objects iterated over ust be printable by std::ofstream
  class Measurer_Collection
  {
  public:
    // Construct given output filename,
    // output scientific notation precision,
    // and output file delimiter
    Measurer_Collection
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
    }
      
    // Destructor cleans up output object
    ~Measurer_Collection()
    { output.close(); }
    
    // Output collection info
    template <typename Subject>
    void operator()(Subject const& subject)
    {
      std::string delim = "";
      for (std::size_t ii = 0; ii < subject.size(); ++ii)
      {
        output << delim << subject.size(ii);
        for (auto it = subject.cbegin(ii); it != subject.cend(ii); ++it)
          output << delimiter << *it;
        delim = delimiter;
      }
      output << "\n";
    }
      
    // Output collection info preceded by scalar value
    template <typename Subject, typename Value>
    void operator()(Subject const& subject, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject);
    }
    
  private:
    std::ofstream output;   // Handle output to file
    std::string delimiter;  // Delimiter between state quantities in each file line
  };
  
  // Store sum of wanted quantity over all particles at each measurement
  // Getter objects in methods must take a particle and output wanted quantity
  // Wanted quantity must be summable by useful::plus_InPlace
  // and printable by useful::print
  // Subject must implement cbegin() and size()
  template <typename Type = double>
  class Measurer_Store_Total
  {
  public:
    // Construct given space to reserve for storing measurements
    Measurer_Store_Total(std::size_t reserve = 0)
    { values.reserve(reserve); }
    
    // Store current total
    template <typename Subject, typename Getter>
    void update(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
        return;
      
      values.push_back(get(*subject.cbegin()));
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend(); ++part_it)
        operation::plus_InPlace(values.back(), get(*part_it));
    }
    
    // Get stored values
    auto const& get()
    { return values; }
    
    // Print totals
    void print
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    {
      std::ofstream output{ filename };
      output << std::setprecision(precision)
             << std::scientific;
      for (auto const& val : values)
        useful::print(output, val, 0, delimiter);
      output << "\n";
      output.close();
    }
    
    // Print each total preceded by each scalar value in measure_points
    template <typename Container>
    void print
    (std::string const& filename, Container const& measure_points,
     int precision = 8, std::string delimiter = "\t")
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < measure_points.size(); ++mm)
      {
        output << measure_points[mm];
        useful::print(output, values[mm], 1, delimiter);
      }
      output << "\n";
      output.close();
    }
    
  private:
    std::vector<Type> values;  // Stored measurements
  };
  
  // Store arithmetic mean of wanted quantity over all particles
  // Getter objects in methods must take a particle and output wanted quantity
  // Wanted quantity must be summable by useful::plus_InPlace,
  // divisible by scalar by useful::div_scalar_InPlace,
  // and printable by useful::print
  // Subject must implement cbegin() and size()
  template <typename Type = double>
  class Measurer_Store_Mean
  {
  public:
    // Construct given space to reserve for storing measurements
    Measurer_Store_Mean(std::size_t reserve = 0)
    { values.reserve(reserve); }
    
    // Store current mean
    template <typename Subject, typename Getter>
    void update(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
        return;
      
      values.push_back(get(*subject.cbegin()));
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend();++part_it)
        operation::plus_InPlace(values.back(), get(*part_it));
      operation::div_scalar_InPlace(values.back(), subject.size());
    }
    
    // Get stored values
    auto const& get()
    { return values; }
    
    // Print means
    void print
    (std::string const& filename, int precision = 8,
     std::string delimiter = "\t")
    {
      std::ofstream output{ filename };
      output << std::setprecision(precision)
             << std::scientific;
      for (auto const& val : values)
      {
        useful::print(output, val, 0, delimiter);
        output << "\n";
      }
      output.close();
    }
    
    // Print each mean preceded by each scalar value in measure_points
    template <typename Container>
    void print
    (std::string const& filename, Container const& measure_points,
     int precision = 8, std::string delimiter = "\t")
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < measure_points.size(); ++mm)
      {
        output << measure_points[mm];
        useful::print(output, values[mm], 1, delimiter);
        output << "\n";
      }
      output.close();
    }
    
  private:
    std::vector<Type> values;
  };
  
  // Store quantity when given values (e.g. of position) are crossed
  // between the new particle state and the old particle state
  // Getter objects in methods must take a particle and output wanted quantity
  // Getter_Crossing objects in methods must take a state and output a position
  // Wanted quantity must be printable by useful::print
  // Subject must implement particles()
  template <typename Type = double>
  class Measurer_Store_Crossing
  {
  public:
    // Construct given crossing values and space to reserve to
    // store measurements for each crossing value
    Measurer_Store_Crossing
    (std::vector<double> crossing_values, std::size_t reserve = 0)
    : crossing_values{ crossing_values }
    , values(crossing_values.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }
    
    // Construct given crossing values and space to reserve to
    // store measurements for each crossing value
    template <typename Container>
    Measurer_Store_Crossing
    (Container const& crossing_values, std::size_t reserve = 0)
    : crossing_values{ get_crossing_values(crossing_values) }
    , values(crossing_values.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }

    // Store new values of each crossing value
    template <typename Subject, typename Getter,
    typename Getter_Crossing = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Crossing const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(crossing_values);
             mm < current_crossed - std::begin(crossing_values); ++mm)
          values[mm].push_back(get(part));
      }
    }

    // Get number of stored values at the mmth crossing value
    std::size_t size(std::size_t mm) const
    { return values[mm].size(); }
    
    // Get stored values
    // rows: crossing value, columns: crossing event
    auto const& get() const
    { return values; }

    // Print crossing value followed by quantity values in each line
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < crossing_values.size(); ++mm)
      {
        output << crossing_values[mm];
        useful::print(output, values[mm], 1, delimiter);
        output << "\n";
      }
      output.close();
    }

    const std::vector<double> crossing_values;  // Store if each of these values is crossed

  private:
    std::vector<std::vector<Type>> values;  // values[ii][jj] is the jjth stored value
                                            // at the iith crossing value
    
    // Convert crossing value container to vector
    template <typename Container>
    auto get_crossing_values(Container const& crossing_values_in)
    {
      std::vector<double> crossing_values;
      for (auto const& val : crossing_values_in)
        crossing_values.push_back(val);
      return crossing_values;
    }
  };
  
  // Store sum of quantity when given values (e.g. of position) are crossed
  // between the new particle state and the old particle state
  // Getter objects in methods must take a particle and output wanted quantity
  // Getter_Crossing objects in methods must take a state and output a position
  // Wanted quantity must be summable by useful::plus_InPlace and
  // printable by useful::print
  // Subject must implement particles()
  template <typename Type = double>
  class Measurer_Store_Crossing_Total
  {
  public:
    // Construct given crossing values
    Measurer_Store_Crossing_Total
    (std::vector<double> crossing_values)
    : crossing_values{ crossing_values }
    , values(crossing_values.size())
    , nr_counts(crossing_values.size())
    {}
    
    // Construct given crossing values
    template <typename Container>
    Measurer_Store_Crossing_Total
    (Container const& crossing_values)
    : crossing_values{ get_crossing_values(crossing_values) }
    , values(crossing_values.size())
    , nr_counts(crossing_values.size())
    {}

    // Update total for each crossing value
    template <typename Subject, typename Getter,
    typename Getter_Crossing = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Crossing const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(crossing_values);
             mm < current_crossed - std::begin(crossing_values); ++mm)
        {
          if (nr_counts[mm] == 0)
            values[mm] = get(part);
          else
            operation::plus_InPlace(values[mm], get(part));
          ++nr_counts[mm];
        }
      }
    }
    
    // Number of crossings for mmth crossing value
    std::size_t counts(std::size_t mm) const
    { return nr_counts[mm]; }
    
    // Get stored totals for each crossing
    auto const& get() const
    { return values; }

    // Print crossing value followed by total value in each line
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < crossing_values.size(); ++mm)
      {
        output << crossing_values[mm];
        useful::print(output, values[mm], 1, delimiter);
        output << "\n";
      }
      output.close();
    }

    const std::vector<double> crossing_values;  // Store if each of these values is crossed

  private:
    std::vector<Type> values;            // Stored totals for each crossing value
    std::vector<std::size_t> nr_counts;  // Number of crossings for each crossing value
    
    // Convert crossing value container to vector
    template <typename Container>
    auto get_crossing_values(Container const& crossing_values_in)
    {
      std::vector<double> crossing_values;
      for (auto const& val : crossing_values_in)
        crossing_values.push_back(val);
      return crossing_values;
    }
  };
      
  // Store quantity when given values (e.g. of position)
  // are crossed for the first time by each particle
  // between the new particle state and the old particle state
  // Getter objects in methods must take a particle and output wanted quantity
  // Getter_Crossing objects in methods must take a state and output a position
  // Wanted quantity must be printable by useful::print
  // Subject must implement particles()
  template <typename Type = double>
  class Measurer_Store_FirstCrossing
  {
  public:
    // Construct given crossing values and space to reserve to
    // store measurements for each crossing value
    Measurer_Store_FirstCrossing
    (std::vector<double> crossing_values, std::size_t reserve = 0)
    : crossing_values{ crossing_values }
    , values(crossing_values.size())
    , particles_crossed(crossing_values.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }
    
    // Construct given crossing values and space to reserve to
    // store measurements for each crossing value
    template <typename Container>
    Measurer_Store_FirstCrossing
    (Container const& crossing_values, std::size_t reserve = 0)
    : crossing_values{ get_crossing_values(crossing_values) }
    , values(crossing_values.size())
    , particles_crossed(crossing_values.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }

    // Store new values for each crossing value
    template <typename Subject, typename Getter,
    typename Getter_Crossing = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Crossing const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(crossing_values);
             mm < current_crossed - std::begin(crossing_values); ++mm)
          if (particles_crossed[mm].insert(part.state_new().tag).second)
            values[mm].push_back(get(part));
      }
    }

    // Get number of particles that have crossed mmth crossing value
    std::size_t size(std::size_t mm) const
    { return values[mm].size(); }
    
    // Get stored values
    // rows: crossing value, columns: crossing event
    auto const& get() const
    { return values; }
    
    // Check if mmth crossing value has been crossed by partth particle
    bool crossed(std::size_t mm, std::size_t part) const
    { return particles_crossed[mm].count(part); };

    // Print crossing value followed by quantity values in each line
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < crossing_values.size(); ++mm)
      {
        output << crossing_values[mm];
        useful::print(output, values[mm], 1, delimiter);
        output << "\n";
      }
      output.close();
    }

    const std::vector<double> crossing_values;  // Store if each of these values is crossed
                                                // for the first time

  private:
    std::vector<std::vector<Type>> values;        // Stored values for each crossing value
    std::vector<std::unordered_set<std::size_t>>
      particles_crossed;                          // Particles that have already crossed
                                                  // each crossing value
    
    // Convert crossing value container to vector
    template <typename Container>
    auto get_crossing_values(Container const& crossing_values_in)
    {
      std::vector<double> crossing_values;
      for (auto const& val : crossing_values_in)
        crossing_values.push_back(val);
      return crossing_values;
    }
  };
      
        
  // Store sum of quantity when given values (e.g. of position)
  // are crossed for the first time by each particle
  // between the new particle state and the old particle state
  // Getter objects in methods must take a particle and output wanted quantity
  // Getter_Crossing objects in methods must take a state and output a position
  // Wanted quantity must be summable by useful::plus_InPlace and
  // printable by useful::print
  // Subject must implement particles()
  template <typename Type = double>
  class Measurer_Store_FirstCrossing_Total
  {
  public:
    // Construct given crossing values
    Measurer_Store_FirstCrossing_Total
    (std::vector<double> crossing_values)
    : crossing_values{ crossing_values }
    , values(crossing_values.size())
    , particles_crossed(crossing_values.size())
    {}
    
    // Construct given crossing values
    template <typename Container>
    Measurer_Store_FirstCrossing_Total
    (Container const& crossing_values)
    : crossing_values{ get_crossing_values(crossing_values) }
    , values(crossing_values.size())
    , particles_crossed(crossing_values.size())
    {}

    // Update total for each crossing value
    template <typename Subject, typename Getter,
    typename Getter_Crossing = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Crossing const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(crossing_values),
          std::end(crossing_values), get_position(part.state_new()));
        
        for (std::size_t mm = last_crossed - std::begin(crossing_values);
             mm < current_crossed - std::begin(crossing_values); ++mm)
          if (particles_crossed[mm].insert(part.state_new().tag).second)
          {
            if (counts(mm) == 1)
              values[mm] = get(part);
            else
              operation::plus_InPlace(values[mm], get(part));
          }
      }
    }
    
    // Get number of particles that have crossed mmth crossing value
    std::size_t counts(std::size_t mm) const
    { return particles_crossed[mm].size(); }
    
    // Get stored totals for each crossing value
    auto const& get() const
    { return values; }
    
    // Check if mmth crossing value has been crossed by partth particle
    bool crossed(std::size_t mm, std::size_t part) const
    { return particles_crossed[mm].count(part); };

    // Print crossing value followed by total value in each line
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t mm = 0; mm < crossing_values.size(); ++mm)
      {
        output << crossing_values[mm];
        useful::print(output, values[mm], 1, delimiter);
        output << "\n";
      }
      output.close();
    }

    const std::vector<double> crossing_values;  // Store if each of these values is crossed
                                                // for the first time

  private:
    std::vector<Type> values;                     // Stored total for each crossing value
    std::vector<std::unordered_set<std::size_t>>
      particles_crossed;                          // Particles that have already crossed
                                                  // each crossing value
    
    // Convert crossing value container to vector
    template <typename Container>
    auto get_crossing_values(Container const& crossing_values_in)
    {
      std::vector<double> crossing_values;
      for (auto const& val : crossing_values_in)
        crossing_values.push_back(val);
      return crossing_values;
    }
  };
      
  // Store change in quantity between returns to a given state quantity
  // ReturnCriterium objects in methods must take a state
  // and return true if a criterium is met and fals otherwise
  // Getter objects in methods must take a state and return the value of the quantity
  // whose change is to be measured between returns
  // Wanted quantity change must be a double
  // Subject must implement particles()
  class Measurer_Store_Return
  {
  public:
    // Construct given subject containing particles,
    // space to reserve for return value measurements, and getter
    // to get return values
    template <typename Subject, typename Getter = ctrw::Get_time>
    Measurer_Store_Return
    (Subject const& subject,
     std::size_t reserve = 0, Getter const& get = {})
    {
      values.reserve(reserve);
      last_left.reserve(subject.size());
      for (auto const& part : subject.particles())
        last_left.push_back(get(part.state_new()));
    }

    // Store return value changes for particles that returned
    template <typename Subject, typename ReturnCriterium,
    typename Getter = ctrw::Get_time>
    void update
    (Subject const& subject, ReturnCriterium const& criterium,
     Getter const& get = {})
    {
      for (auto const& part : subject.particles())
      {
        auto const& state_old = part.state_old();
        auto const& state_new = part.state_new();
        bool was_inside = criterium(state_old);
        bool is_inside = criterium(state_new);
        if (is_inside && !was_inside)
          values.push_back(get(state_new)-last_left[state_new.tag]);
        else if (!is_inside && was_inside)
          last_left[state_new.tag] = get(state_new);
      }
    }
    
    // Get stored number of returns
    std::size_t counts() const
    { return values.size(); }
    
    // Get stored return value changes
    auto const& get() const
    { return values; }
    
    // Output stored return value changes
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(precision)
             << std::scientific;
      useful::print(output, values);
      output << "\n";
      output.close();
    }
    
  private:
    std::vector<double> values;     // Stored return value changes
    std::vector<double> last_left;  // Value of quantity at last return
                                    // for each particle
  };
}

#endif /* Measurer_CTRW_h */
