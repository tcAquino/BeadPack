//
//  ParticleCollections.h
//  BeadPack
//
//  Created by Tomás Aquino on 18/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef ParticleCollections_h
#define ParticleCollections_h

#include <iterator>
#include <list>
#include "StateGetter.h"
#include "general/Operations.h"
#include "general/useful.h"

namespace ctrw
{
  // Note: For a closed strip (loop), repeat the 1st particle at the end
  // Warning: Requires particle tags to match position in CTRW particle container.
  class Strip
  {
  private:
    using Container = std::list<std::size_t>;
    using iterator = Container::iterator;
    
  public:
    using const_iterator = Container::const_iterator;
    
    template <typename Container = std::vector<std::size_t>>
    Strip
    (Container const& particle_tags_initial, std::size_t max_particles, double max_segment_length)
    : max_particles{ max_particles }
    , max_segment_length{ max_segment_length }
    {
      for (auto const& tag : particle_tags_initial)
        particle_tags.push_back(tag);
    }
    
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void resize(CTRW& ctrw, Getter_position get_position, Adjust adjust, StateMaker state_maker)
    {
      for (auto it = std::next(begin()); it != end();)
      {
        if (full())
          return;
        if (segment_length_sq(it, ctrw, get_position) > max_segment_length*max_segment_length)
        {
          insert_midpoint(it, ctrw, get_position, adjust, state_maker);
          std::advance(it, -1);
        }
        else
        {
          ++it;
        }
      }
    }
    
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void refine_with_old_state(CTRW& ctrw, Getter_position get_position, Adjust adjust, StateMaker state_maker)
    {
      for (auto it = std::next(begin()); it != end(); ++it)
      {
        if (full())
          return;
        if (segment_length_sq(it, ctrw, get_position) > max_segment_length*max_segment_length)
          insert_previous_state(it, ctrw, get_position, adjust, state_maker);
      }
    }
    
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void refine(std::size_t goal_size, CTRW& ctrw, Getter_position get_position, Adjust adjust, StateMaker state_maker)
    {
      while (1)
        for (auto it = std::next(begin()); it != end(); ++it)
        {
          if (full() || size() >= goal_size)
            return;
          insert_midpoint(it, ctrw, get_position, adjust, state_maker);
        }
    }
    
    template <typename CTRW, typename Getter_position>
    double length(CTRW const& ctrw, Getter_position get_position) const
    {
      double len = 0.;
      for (auto it = std::next(cbegin()); it != cend(); ++it)
        len += segment_length(it, ctrw, get_position);
      
      return len;
    }
    
    std::size_t size() const
    { return particle_tags.size(); }
    
    std::size_t maximum_size() const
    { return max_particles; }
    
    bool full() const
    { return size() >= maximum_size(); }
    
    const_iterator cbegin() const
    { return particle_tags.cbegin(); }
    
    const_iterator cend() const
    { return particle_tags.cend(); }
    
  private:
    std::size_t max_particles;
    double max_segment_length;
    Container particle_tags;
    
    template <typename CTRW, typename Getter_position>
    double segment_length(iterator it, CTRW const& ctrw, Getter_position get_position) const
    {
      return
      operation::abs(operation::minus(get_position(ctrw.particles(*it)),
                                      get_position(ctrw.particles(*std::prev(it)))));
    }
    
    template <typename CTRW, typename Getter_position>
    double segment_length_sq(iterator it, CTRW const& ctrw, Getter_position get_position) const
    {
      return
      operation::abs_sq(operation::minus(get_position(ctrw.particles(*it)),
                                         get_position(ctrw.particles(*std::prev(it)))));
    }
    
    template <typename CTRW, typename Getter_position>
    auto segment_midpoint(iterator it, CTRW const& ctrw, Getter_position get_position) const
    {
      auto midpoint = operation::plus(get_position(ctrw.particles(*it)),
                                      get_position(ctrw.particles(*std::prev(it))));
      operation::times_scalar_InPlace(0.5, midpoint);
      
      return midpoint;
    }
    
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void insert_midpoint
    (iterator it, CTRW& ctrw, Getter_position get_position,
     Adjust adjust, StateMaker state_maker)
    {
      auto state{ state_maker() };
      state.position = segment_midpoint(it, ctrw, get_position);
      state.tag = ctrw.size();
      adjust(state);
      
      ctrw.push_back(Particle{ state });
      particle_tags.insert(it, state.tag);
    }
    
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void insert_previous_state
    (iterator it, CTRW& ctrw, Getter_position get_position,
     Adjust adjust, StateMaker state_maker)
    {
      particle_tags.insert(it, ctrw.size());
      ctrw.push_back(Particle{ ctrw.particles(*it).State_old() });
    }
    
    iterator begin()
    { return particle_tags.begin(); }
    
    iterator end()
    { return particle_tags.end(); }
  };
  
  template
  <typename CTRW, typename Getter_position = ctrw::Get_new_from_particle<Get_position>,
  typename Adjust = useful::DoNothing,
  typename StateMaker = useful::Maker<typename CTRW::State>>
  struct StripHandler
  {
    CTRW& ctrw;
    Getter_position get_position;
    Adjust adjust;
    StateMaker state_maker;
    
    std::vector<Strip> strips;
    
    StripHandler
    (CTRW& ctrw,
     Getter_position get_position = {},
     Adjust adjust = {}, StateMaker state_maker = {})
    : ctrw{ ctrw }
    , get_position{ get_position }
    , adjust{ adjust }
    , state_maker{ state_maker }
    {}
    
    void reserve(std::size_t nr_strips)
    { strips.reserve(nr_strips); }
    
    std::size_t size() const
    { return strips.size(); }
    
    template <typename Container = std::vector<std::size_t>>
    void push_back(Container const& particle_tags_initial, std::size_t max_particles, double max_segment_length)
    {
      strips.emplace_back(particle_tags_initial, max_particles, max_segment_length);
    }
    
    void refine(std::size_t goal_size)
    {
      for (auto& strip : strips)
        strip.refine(goal_size, ctrw, get_position, adjust, state_maker);
    }
    
    void refine(std::size_t ss, std::size_t goal_size)
    {
      strips[ss].refine(goal_size, ctrw, get_position, adjust, state_maker);
    }
    
    void resize(std::size_t ss)
    {
      strips[ss].resize(ctrw, get_position, adjust, state_maker);
    }
    
    void resize()
    {
      for (auto& strip : strips)
        strip.resize(ctrw, get_position, adjust, state_maker);
    }
    
    void refine_with_old_state(std::size_t ss)
    {
      strips[ss].refine_with_old_state(ctrw, get_position, adjust, state_maker);
    }
    
    void refine_with_old_state()
    {
      for (auto& strip : strips)
        strip.refine_with_old_state(ctrw, get_position, adjust, state_maker);
    }
    
    std::size_t size(std::size_t ss) const
    { return strips[ss].size(); }
    
    double length(std::size_t ss) const
    { return strips[ss].length(ctrw, get_position); }
    
    std::size_t maximum_size(std::size_t ss) const
    { return strips[ss].maximum_size; }
    
    bool full(std::size_t ss) const
    { return strips[ss].full(); }
    
    bool full() const
    {
      for (auto strip : strips)
        if (!strip.full())
          return 0;
      return 1;
    }
    
    auto cbegin(std::size_t ss) const
    { return strips[ss].cbegin(); }
    
    auto cend(std::size_t ss) const
    { return strips[ss].cend(); }
  };
}







#endif /* ParticleCollections_h */
