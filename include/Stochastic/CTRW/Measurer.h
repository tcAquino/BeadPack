//
//  Measurer.h
//  CTRW_velocity
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
  class Measurer_State
  {
  public:
    struct New{};
    struct Old{};
    struct Both{};
    
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
    
    ~Measurer_State()
    { output.close(); }
      
    template <typename CTRW, typename OutputState>
    void operator()
    (CTRW const& ctrw, OutputState output_state, New)
    {
      bool delim = 0;
      for (auto const& part : ctrw.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_new(), output, delimiter);
        delim = 1;
      }
      output << "\n";
    }
      
    template <typename CTRW, typename OutputState>
    void operator()
    (CTRW const& ctrw, OutputState output_state, Old)
    {
      bool delim = 0;
      for (auto const& part : ctrw.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_old(), output, delimiter);
        delim = 1;
      }
      output << "\n";
    }
    
    template <typename CTRW, typename OutputState>
    void operator()
    (CTRW const& ctrw, OutputState output_state, Both)
    {
      bool delim = 0;
      for (auto const& part : ctrw.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_new(), output, delimiter);
        output << delimiter;
        output_state(part.state_old(), output, delimiter);
      }
      output << "\n";
    }
      
    template <typename CTRW, typename OutputState,
    typename Value,typename Which>
    void operator()
    (CTRW const& ctrw, OutputState output_state, Value tag,
     Which which)
    {
      output << tag << delimiter;
      (*this)(ctrw, output_state, which);
    }
    
  private:
    std::ofstream output;
    std::string delimiter;
  };
      
  class Measurer_Collection
  {
  public:
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
      
    ~Measurer_Collection()
    { output.close(); }
    
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
      
    template <typename Subject, typename Value>
    void operator()(Subject const& subject, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject);
    }
    
  private:
    std::ofstream output;
    const std::string delimiter;
  };
  
  class Measurer_Particle
  {
  public:
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
    
    ~Measurer_Particle()
    { output.close(); }
    
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
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
  private:
    std::ofstream output;
    const std::string delimiter;
  };
      
  class Measurer_Total
  {
  public:
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
    
    ~Measurer_Total()
    { output.close(); }
    
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
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
    private:
      std::ofstream output;
      const std::string delimiter;
  };
      
  class Measurer_Mean
  {
  public:
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
    
    ~Measurer_Mean()
    { output.close(); }
    
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
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
    private:
      std::ofstream output;
      const std::string delimiter;
  };
      
  template <typename Type = double>
  class Measurer_Store_Total
  {
  public:
    Measurer_Store_Total(std::size_t reserve = 0)
    { values.reserve(reserve); }
    
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
      {
        values.push_back();
        return;
      }
      auto val = get(*subject.cbegin());
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend(); ++part_it)
        operation::plus_InPlace(val, get(*part_it));
    }
    
    auto const& get()
    { return values; }
    
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
    
    template <typename Cont>
    void print
    (std::string const& filename, Cont const& measure_points,
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
    std::vector<Type> values;
  };
      
  template <typename Type = double>
  class Measurer_Store_Mean
  {
  public:
    Measurer_Store_Mean(std::size_t reserve = 0)
    { values.reserve(reserve); }
    
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      if (subject.size() == 0)
      {
        values.push_back();
        return;
      }
      
      auto val = get(*subject.cbegin());
      for (auto part_it = std::next(subject.cbegin());
           part_it != subject.cend();++part_it)
        operation::plus_InPlace(val, get(*part_it));
      operation::div_scalar_InPlace(val, subject.size());
    }
    
    auto const& get()
    { return values; }
    
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
    
    template <typename Cont>
    void print
    (std::string const& filename, Cont const& measure_points,
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
  
  template <typename Type = double>
  class Measurer_Store_Crossing
  {
  public:
    Measurer_Store_Crossing
    (std::vector<double> measure_points, std::size_t reserve = 0)
    : measure_points{ measure_points }
    , values(measure_points.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }
    
    template <typename Container>
    Measurer_Store_Crossing
    (Container const& measure_points, std::size_t reserve = 0)
    : measure_points{ get_measure_points(measure_points) }
    , values(measure_points.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }

    template <typename Subject, typename Getter,
    typename Getter_Position = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Position const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(measure_points);
             mm < current_crossed - std::begin(measure_points); ++mm)
          values[mm].push_back(get(part));
      }
    }

    std::size_t size(std::size_t mm) const
    { return values[mm].size(); }
    
    auto const& get() const
    { return values; }

    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
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

    const std::vector<double> measure_points;

  private:
    std::vector<std::vector<Type>> values;
    
    template <typename Container>
    auto get_measure_points(Container const& measure_points_in)
    {
      std::vector<double> measure_points;
      for (auto const& val : measure_points_in)
        measure_points.push_back(val);
      return measure_points;
    }
  };
      
  template <typename Type = double>
  class Measurer_Store_Crossing_Total
  {
  public:
    Measurer_Store_Crossing_Total
    (std::vector<double> measure_points)
    : measure_points{ measure_points }
    , values(measure_points.size())
    , nr_counts(measure_points.size())
    {}
    
    template <typename Container>
    Measurer_Store_Crossing_Total
    (Container const& measure_points)
    : measure_points{ get_measure_points(measure_points) }
    , values(measure_points.size())
    , nr_counts(measure_points.size())
    {}

    template <typename Subject, typename Getter,
    typename Getter_Position = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Position const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(measure_points);
             mm < current_crossed - std::begin(measure_points); ++mm)
        {
          if (nr_counts[mm] == 0)
            values[mm] = get(part);
          else
            operation::plus_InPlace(values[mm], get(part));
          ++nr_counts[mm];
        }
      }
    }
    
    std::size_t counts(std::size_t mm) const
    { return nr_counts[mm]; }
    
    auto const& get() const
    { return values; }

    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
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

    const std::vector<double> measure_points;

  private:
    std::vector<Type> values;
    std::vector<std::size_t> nr_counts;
    
    template <typename Container>
    auto get_measure_points(Container const& measure_points_in)
    {
      std::vector<double> measure_points;
      for (auto const& val : measure_points_in)
        measure_points.push_back(val);
      return measure_points;
    }
  };
      
  template <typename Type = double>
  class Measurer_Store_FirstCrossing
  {
  public:
    Measurer_Store_FirstCrossing
    (std::vector<double> measure_points, std::size_t reserve = 0)
    : measure_points{ measure_points }
    , values(measure_points.size())
    , particles_crossed(measure_points.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }
    
    template <typename Container>
    Measurer_Store_FirstCrossing
    (Container const& measure_points, std::size_t reserve = 0)
    : measure_points{ get_measure_points(measure_points) }
    , values(measure_points.size())
    , particles_crossed(measure_points.size())
    {
      for (auto& val : values)
        val.reserve(reserve);
    }

    template <typename Subject, typename Getter,
    typename Getter_Position = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Position const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_new()));
        for (std::size_t mm = last_crossed - std::begin(measure_points);
             mm < current_crossed - std::begin(measure_points); ++mm)
          if (particles_crossed[mm].insert(part.state_new().tag).second)
            values[mm].push_back(get(part));
      }
    }

    std::size_t size(std::size_t mm) const
    { return values[mm].size(); }
    
    auto const& get() const
    { return values; }
    
    bool crossed(std::size_t mm, std::size_t part) const
    { return particles_crossed[mm].count(part); };

    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
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

    const std::vector<double> measure_points;

  private:
    std::vector<std::vector<Type>> values;
    std::vector<std::unordered_set<std::size_t>> particles_crossed;
    
    template <typename Container>
    auto get_measure_points(Container const& measure_points_in)
    {
      std::vector<double> measure_points;
      for (auto const& val : measure_points_in)
        measure_points.push_back(val);
      return measure_points;
    }
  };
      
  template <typename Type = double>
  class Measurer_Store_FirstCrossing_Total
  {
  public:
    Measurer_Store_FirstCrossing_Total
    (std::vector<double> measure_points)
    : measure_points{ measure_points }
    , values(measure_points.size())
    , particles_crossed(measure_points.size())
    {}
    
    template <typename Container>
    Measurer_Store_FirstCrossing_Total
    (Container const& measure_points)
    : measure_points{ get_measure_points(measure_points) }
    , values(measure_points.size())
    , particles_crossed(measure_points.size())
    {}

    template <typename Subject, typename Getter,
    typename Getter_Position = ctrw::Get_position_component<0>>
    void update
    (Subject const& subject, Getter const& get,
     Getter_Position const& get_position = {})
    {
      for (auto const& part : subject.particles())
      {
        auto last_crossed = std::upper_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_old()));
        auto current_crossed = std::lower_bound(std::begin(measure_points),
          std::end(measure_points), get_position(part.state_new()));
        
        for (std::size_t mm = last_crossed - std::begin(measure_points);
             mm < current_crossed - std::begin(measure_points); ++mm)
          if (particles_crossed[mm].insert(part.state_new().tag).second)
          {
            if (counts(mm) == 1)
              values[mm] = get(part);
            else
              operation::plus_InPlace(values[mm], get(part));
          }
      }
    }
    
    std::size_t counts(std::size_t mm) const
    { return particles_crossed[mm].size(); }
    
    auto const& get() const
    { return values; }
    
    bool crossed(std::size_t mm, std::size_t part) const
    { return particles_crossed[mm].count(part); };

    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
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

    const std::vector<double> measure_points;

  private:
    std::vector<Type> values;
    std::vector<std::unordered_set<std::size_t>> particles_crossed;
    
    template <typename Container>
    auto get_measure_points(Container const& measure_points_in)
    {
      std::vector<double> measure_points;
      for (auto const& val : measure_points_in)
        measure_points.push_back(val);
      return measure_points;
    }
  };
      
  class Measurer_Store_Return
  {
  public:
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
    
    std::size_t counts() const
    { return values.size(); }
    
    auto const& get() const
    { return values; }
    
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
    std::vector<double> values;
    std::vector<double> last_left;
  };
}

#endif /* Measurer_CTRW_h */
