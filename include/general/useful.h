//
// useful.h
//
// Created on: Mar 15, 2011
// Author: tomas
//

// Miscelaneous collection of useful objects and algorithms

#ifndef USEFUL_H_
#define USEFUL_H_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace useful
{
  
  // Check if container contains val
  // Warning: container must be sorted
  template <typename T, typename U>
  bool contains(const std::vector<T>& container, U const& val)
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          [](T const& elem, U const& val){ return elem < val; });
    
    return it != container.end() && *it == val;
  }
  
  // Split string
  // Adapted from Beder Acosta Borges's answer here:
  // https://stackoverflow.com/questions/14265581/
  // parse-split-a-string-in-c-using-string-delimiter-standard-c
  bool endsWith(const std::string& string, const std::string& suffix)
  {
    return string.size() >= suffix.size() &&
      string.substr(string.size() - suffix.size()) == suffix;
  }
  
  std::vector<std::string> split
  (std::string const& string, std::string const& delimiter = " ",
   bool empty_entries = false)
  {
    std::vector<std::string> tokens;

    for (std::size_t start = 0, end; start < string.length();
         start = end+delimiter.length())
    {
      std::size_t position = string.find(delimiter, start);
      end = position != std::string::npos? position : string.length();

      std::string token = string.substr(start, end-start);
      if (empty_entries || !token.empty())
        tokens.push_back(token);
    }

    if (empty_entries &&
        (string.empty() || endsWith(string, delimiter)))
      tokens.push_back("");

    return tokens;
  }
  
  auto parse_error
  (std::string const& filename, std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line +
      "\nin file " + filename };
  }
  
  auto parse_error_file(std::string const& filename)
  {
    return std::runtime_error{
      "Could not parse file " + filename };
  }
  
  auto parse_error_line(std::string const& line)
  {
    return std::runtime_error{
      "Could not parse line\n" + line };
  }
  
  auto open_read_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for reading" };
  }
  
  auto open_write_error(std::string const& filename)
  {
    return std::runtime_error{
      "Could not open file " + filename + " for writing" };
  }
  
  auto bad_file_contents(std::string const& filename)
  {
    return std::runtime_error{
      "Innapropriate contents in file " + filename };
  }
  
  auto bad_eof
  (std::string const& filename, std::string const& string)
  {
    return std::runtime_error{
    "Reached end of " +
      filename + " before " + string + " was found" };
  }
  
  auto bad_parameters()
  {
    return std::invalid_argument{ "Inappropriate parameters" };
  }
  
  //Load 1-column file into vectors of doubles
  auto load_1
  (std::string const& filename, std::size_t nr_estimate = 0,
   std::size_t header_lines = 0,
   std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    Container values;
    values.reserve(nr_estimate);
    
    std::ifstream file(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != 1)
        throw parse_error(filename, line);
      values.push_back(std::stod(split_line[0]));
    }
    file.close();
    
    return values;
  }
  
  //Load 2-column file into pair of vectors of doubles
  auto load_2
 (std::string const& filename, std::size_t nr_estimate = 0,
  std::size_t header_lines = 0, std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    std::pair<Container, Container> values;
    values.first.reserve(nr_estimate);
    values.second.reserve(nr_estimate);
    
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != 2)
        throw parse_error(filename, line);
      
      values.first.push_back(std::stod(split_line[0]));
      values.second.push_back(std::stod(split_line[1]));
    }
    file.close();
    
    return values;
  }
  
  //Load file into vector of vectors of doubles
  auto load(std::string const& filename, std::size_t nr_columns,
            std::size_t nr_estimate = 0,
              std::size_t header_lines = 0,
              std::string const& delim = " ")
  {
    using Value = double;
    using Container = std::vector<Value>;
    std::vector<Container> values(nr_columns);
    for (auto& val : values)
      val.reserve(nr_estimate);
    
    std::ifstream file(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line = split(line, delim);
      if (split_line.size() != nr_columns)
        throw parse_error(filename, line);
      
      for (std::size_t cc = 0; cc < nr_columns; ++cc)
        values[cc].push_back(std::stod(split_line[cc]));
    }
    file.close();
    
    return values;
  }
  
  // From Anton Dyachenko's answer here:
  // https://stackoverflow.com/questions/6534041/how-to-check-whether-operator-exists
  template <class T, class R, class ... Args>
  std::is_convertible<std::invoke_result_t<T, Args...>, R> is_invokable_test(int);
  template <class T, class R, class ... Args>
  std::false_type is_invokable_test(...);
  template <class T, class R, class ... Args>
  using is_invokable = decltype(is_invokable_test<T, R, Args...>(0));
  template <class T, class R, class ... Args>
  constexpr auto is_invokable_v = is_invokable<T, R, Args...>::value;
  
  template <class L, class R = L>
  using has_equality = is_invokable<std::equal_to<>, bool, L, R>;
  template <class L, class R = L>
  constexpr auto has_equality_v = has_equality<L, R>::value;
  
  template <class L, class R = L>
  using has_plus = is_invokable<std::plus<>, bool, L, R>;
  template <class L, class R = L>
  constexpr auto has_plus_v = has_plus<L, R>::value;
  
  template <class L, class R = L>
  using has_minus = is_invokable<std::minus<>, bool, L, R>;
  template <class L, class R = L>
  constexpr auto has_minus_v = has_minus<L, R>::value;
  
  template <class L, class R = L>
  using has_multiplies = is_invokable<std::multiplies<>, bool, L, R>;
  template <class L, class R = L>
  constexpr auto has_multiplies_v = has_multiplies<L, R>::value;
  
  template <class L, class R = L>
  using has_divides = is_invokable<std::divides<>, bool, L, R>;
  template <class L, class R = L>
  constexpr auto has_divides_v = has_divides<L, R>::value;
  
  // From Passer By's answer here:
  // https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain
  template <typename = void, typename... Args>
  struct can_call_sqrt : std::false_type {};
  template <typename... Args>
  struct can_call_sqrt<
  std::void_t<decltype(std::sqrt(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  template <typename... Args>
  inline constexpr bool can_call_sqrt_v = can_call_sqrt<void, Args...>::value;
  
  template <typename = void, typename... Args>
  struct can_call_abs : std::false_type {};
  template <typename... Args>
  struct can_call_abs<
  std::void_t<decltype(std::abs(std::declval<Args>()...))>, Args...>
  : std::true_type {};
  template <typename... Args>
  inline constexpr bool can_call_abs_v = can_call_abs<void, Args...>::value;
  
  template <typename T, typename U, typename Comp_less, typename Comp_eq>
  bool contains
  (const std::vector<T>& container, U const& val,
   Comp_less comp_less, Comp_eq comp_eq)
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          comp_less);
    
    return it != container.end() && comp_eq(*it, val);
  }

  // Sign of val
  template <typename T> int sgn(T val)
  { return (T(0) < val) - (val < T(0)); }

  template <typename Object_Type, typename Return_Type = Object_Type>
  struct StoreConst
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    const Object_Type obj;
  };

  template <typename Object_Type, typename Return_Type = Object_Type>
  struct Store
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    Object_Type obj;
  };

  struct Empty
  {
    template <typename ...Args>
    Empty(Args...){}
    Empty(){}
  };
  
  template <typename Stream, typename Container>
  void print
  (Stream& stream, Container const& container,
   bool delimit_first = 0, std::string delimiter = "\t")
  {
    // TODO: Choose this specialization when stream << container exists
    if constexpr (std::is_pod<Container>::value)
    {
      if (delimit_first)
        stream << delimiter;
      stream << container;
    }
    else
    {
      std::string delim = delimit_first ? delimiter : "";
      for (auto const& val : container)
      {
        stream << delim << val;
        delim = delimiter;
      }
    }
  }
  
  auto read
  (std::string const& filename)
  {
    std::ifstream file{ filename };
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::vector<double> vals;
    double val;
    while(file >> val)
      vals.push_back(val);
    file.close();

    return vals;
  }
  
  struct DoNothing
  {
    template <typename ...Args>
    void operator()(Args...) const
    {}
    
  };

  template <typename ...Args>
  struct DoFalse
  { bool operator()(Args...) const { return false; } };

  template <typename Object, typename Params>
  Object Create(Params const& parameters = useful::Empty())
  {
    Object tracker;
    tracker.Initialize(parameters);

    return tracker;
  }
  
  template <typename Object>
  class Maker
  {
    template <typename ...Args>
    Object operator()(Args... args)
    { return Object{ args... }; }
  };

  template <typename T>
  struct Creator
  {
    typedef T value_type;
    typedef T* pointer;

    pointer operator() (void) const
    { return (new T ()); }

    pointer operator() (double param) const
    { return (new T(param)); }
  };

  template <typename Object_Type>
  struct Forward
  {
    Object_Type operator()(Object_Type const& object)
    { return object; }
  };

  template <typename Object_Type>
  struct Forward_ref
  {
    Object_Type const& operator()(Object_Type const& object)
    { return object; }
  };

  struct Bin
  {
    double left_edge;
    double right_edge;

    double mean() const
    { return (left_edge + right_edge)/2.; }

    double width() const
    { return right_edge - left_edge; }
  };
  
  template
  <typename Value, typename Container_boundaries, typename Container>
  void bin(Value const& value, Container_boundaries const& boundaries, Container& hist)
  {
    if (value < boundaries[0])
    {
      ++hist[0];
      return;
    }
    if (value >= boundaries[boundaries.size()-1])
    {
      ++hist[boundaries.size()-2];
      return;
    }
    
    auto upper_edge = std::upper_bound(std::begin(boundaries),
                                       std::end(boundaries),
                                       value);
    std::size_t idx = upper_edge - std::begin(boundaries);
    ++hist[idx-1];
  }
  
  template
  <typename Value, typename Container_boundaries, typename Container, typename Compare>
  void bin(Value const& value, Container_boundaries const& boundaries, Container const& hist, Compare comp)
  {
    if (comp(value, boundaries[0]))
    {
      ++hist[0];
      return;
    }
    if (!comp(value, boundaries[boundaries.size()-1]))
    {
      ++hist[boundaries.size()-2];
      return;
    }
    
    auto upper_edge = std::upper_bound(std::begin(boundaries),
                                       std::end(boundaries),
                                       value);
    std::size_t idx = upper_edge - std::begin(boundaries);
    ++hist[idx-1];
  }

  template <typename TT> struct Selector_t{};
  template <typename TT, TT(val)> struct Selector{};

  template <typename T>
  T& ensure_ref(T const* obj)
  { return *obj; }

  template <typename T>
  T& ensure_ref(T const& obj)
  { return obj; }

  // Simple hash for a container by combining element hashes
  template <typename Container>
  struct hash_container
  {
    std::size_t operator()(Container const& container) const
    {
      using value_type = typename Container::value_type;
      std::size_t hash_val = std::hash<value_type>()(container[ 0 ]);
      for (size_t ii = 1; ii < container.size() - 1; ++ii)
      {
        hash_val = (hash_val ^ (std::hash<value_type>()(container[ii]) << 1)) >> 1;
      }
      hash_val = hash_val ^ (std::hash<value_type>()(container[container.size()-1]) << 1);

      return hash_val;
    }
  };

  template <typename T>
  void hash_combine(std::size_t & seed, T const& v)
  {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  template <typename S, typename T>
  struct hash_pair
  {
    std::size_t operator()(const std::pair< S, T > & v) const
    {
      std::size_t seed = 0;
      hash_combine(seed, v.first);
      hash_combine(seed, v.second);
      return seed;
    }
  };

  // Replace NaNs
  template < typename T >
  void deNaN(T& to_replace, T replace_with = T())
  {
    if (to_replace != to_replace) to_replace = replace_with;
  }

  // Replace values below tolerance
  template <typename T>
  void chop(T& to_chop, T tolerance, T replace_with = T())
  {
    if (std::abs(to_chop - replace_with) < tolerance) to_chop = replace_with;
  }

  // Copy end element of container to position, then erase end element
  template <typename Container>
  void swap_erase(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    container.pop_back();
  }

  // Copy pointer to end element of vector to position, then delete end element
  template < typename Container >
  void swap_delete(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    delete container.back();
    container.pop_back();
  }
  
  std::size_t countlines(FILE *fin)
  {
    std::size_t lines = 0;
    int maxlength = 255;
    char buffer[ maxlength + 1 ];

    while (fgets(buffer , maxlength + 1 , fin) != NULL)
    {
      if (buffer[std::strlen(buffer) - 1] != '\n')
        throw "Line too long.";
      ++lines;
    }
    if (!feof(fin))
      throw "Could not reach end of file.";

    rewind(fin);

    return lines;
  }

  //	To check if a class has a method
  //	From here: https://stackoverflow.com/questions/29772601/why-is-sfinae-causing-failure-when-there-are-two-functions-with-different-signat
  //	bundle of types
  template <class...> struct types{using type=types;};
  template <class...> struct voider{using type=void;};
  //	Some dists still do not to have void_t
  template <class...Ts> using void_t=typename voider<Ts...>::type;
  //	hide the SFINAE stuff in a details namespace:
  namespace details
  {
    template <template <class...> class Z, class types, class=void>
    struct has_method : std::false_type {};
    template <template <class...> class Z, class...Ts>
    struct has_method<Z,types<Ts...>,void_t<Z<Ts...>>>:std::true_type{};
  }
  // has_method<template, types...> is true iff template <types...> is valid
  template <template <class...> class Z, class...Ts>
  using has_method=details::has_method<Z,types<Ts...>>;

  // Implemented here because some versions of gcc have trouble with the standard one
  template <typename T>
  bool isnan(T const& val)
  { return val != val; }

  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple const& tuple, F f, std::index_sequence< Indices... >)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(tuple)), void(), int{})... };
  }

  template <typename Tuple, typename F>
  void for_each(Tuple const& tuple, F f)
  {
    constexpr std::size_t N = std::tuple_size< std::remove_reference_t< Tuple > >::value;
    for_each_impl(tuple, f, std::make_index_sequence< N >{});
  }

  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple&& tuple, F&& f, std::index_sequence< Indices... >)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get< Indices >(std::forward< Tuple >(tuple))), void(), int{})... };
  }

  template <typename Tuple, typename F>
  void for_each(Tuple&& tuple, F&& f)
  {
    constexpr std::size_t N = std::tuple_size< std::remove_reference_t< Tuple > >::value;
    for_each_impl(std::forward<Tuple>(tuple), std::forward<F>(f),  std::make_index_sequence<N>{});
  }

  // Indices for template metamagic
  template <std::size_t... Indices>
  struct indices
  { using next = indices< Indices..., sizeof...(Indices) >; };

  template <std::size_t size>
  struct build_indices
  { using type = typename build_indices< size - 1 >::type::next; };

  template <>
  struct build_indices< 0 >
  { using type = indices<>; };
}

#endif /* USEFUL_H_ */
