//
//  Boundary.h
//  Boundary
//
//  Created by Tomas Aquino on 8/2/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef Boundary_h
#define Boundary_h

#include <cmath>
#include <numeric>
#include <utility>
#include <vector>
#include "general/Constants.h"
#include "Geometry/Shape.h"
#include "Grid/Grid.h"


// A boundary class must implement the following basic functionality
// class Example_Boundary
// {
// public:
//
//   // Check if position is out of bounds
//   template <typename Position>
//   bool outOfBounds(Position const& position) const
//   {
//     // Return true if boundary condition is to be applied, false otherwise
//   }
//
//   // Enforce boundary condition
//   template <typename State>
//   bool operator()(State& state, State const& state_old) const
//   {
//     // Apply the boundary condition given current and previous state
//     // Return true if anything was done, false otherwise
//   }
// };

namespace boundary
{
  // Return displacement to be added to 1d position to place it within
  // periodic boundaries
	double periodic
 (double position, std::pair<double, double> const& boundaries)
	{
		double box_size = boundaries.second - boundaries.first;
		double pos = position - boundaries.first;
		return -std::floor(pos/box_size)*box_size;
	}
  
  // Return pair of displacement to be added to 1d position
  // and the associated signed number of domain displacements out
  // to place it within periodic boundaries
  std::pair<double, int> periodic_with_outside_info
  (double position, std::pair<double, double> const& boundaries)
  {
    double box_size = boundaries.second - boundaries.first;
    double pos = position - boundaries.first;
    int outside = int(std::floor(pos/box_size));
    return { -outside*box_size, outside };
  }

  // Return displacement to be added to 1d position
  // to reflect it off boundaries
	double reflecting(double position, std::pair<double, double> const& boundaries)
	{
		double box_size = boundaries.second - boundaries.first;
		double pos = position - boundaries.first;
		double nr_boxes_jumped = std::trunc(pos / box_size);
		if (int(nr_boxes_jumped) % 2 == 0)
			return -pos + std::abs(pos - nr_boxes_jumped * box_size);
		else
			return -pos + box_size - std::abs(pos - nr_boxes_jumped * box_size);
	}

  // Return whether 1d position is outside boundaries
  bool outOfBounds_box(double position, std::pair<double,double> const& boundaries)
  {
    return position < boundaries.first || position > boundaries.second;
  }

  // Return whether position is outside boundaries along each dimension
  template <typename Position, typename Boundaries>
  bool outOfBounds_box(Position const& position, Boundaries const& boundaries)
  {
    position[0];
    for (std::size_t dd=0; dd<position.size(); ++dd)
      if (outOfBounds_box(position[dd], boundaries[dd]))
        return true;

    return false;
  }

  // Boundary with no effect
	class DoNothing
	{
	public:
    // Check if position is out of bounds (it never is)
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return false; }

    // Enforce boundary condition (do nothing)
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    { return false; }
	};

  // Periodic boundary in 1d
  // State must define: position (scalar type)
	class Periodic_1d
	{
	public:
		std::pair<double, double> boundaries;  // Lower and upper boundaries

    // Construct given lower and upper boundaries
		Periodic_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

    // Check if position is out of bounds
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return OutOfBounds_Box(position, boundaries); }

    // Enforce boundary condition
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += periodic(state.position, boundaries);
      return true;
    }
	};

  // Periodic boundaries along each dimension
  // State must define: position (vector type)
	class Periodic
	{
	public:
    // Lower and upper boundaries along each dimension
		const std::vector<std::pair<double, double>> boundaries;

    // Construct given lower and upper boundaries along each dimension
		Periodic(std::vector<std::pair<double, double>> boundaries)
		: boundaries(boundaries)
		{}

    // Check if position is out of bounds
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    // Enforce boundary condition
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
        state.position[dd] += periodic(state.position[dd], boundaries[dd]);
      return true;
		}
	};
  
  // Periodic boundaries along each dimension
  // with information about where position would be outside domain
  // State must define: position (vector type)
  //                    periodicity (std::vector<int>) counting how many domains
  //                                                   have been traveled along each dimension
  class Periodic_WithOutsideInfo
  {
  public:
    const std::vector<std::pair<double, double>>
      boundaries;  // Lower and upper boundaries along each dimension

    // Construct given lower and upper boundaries
    Periodic_WithOutsideInfo
    (std::vector<std::pair<double, double>> boundaries)
    : boundaries(boundaries)
    {}

    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position, boundaries); }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
      {
        auto change_outside =
          periodic_with_outside_info(state.position[dd], boundaries[dd]);
        state.position[dd] += change_outside.first;
        state.periodicity[dd] += change_outside.second;
      }
      return true;
    }
  };

  // Reflecting boundary in one dimension
  // State must define: position (scalar type)
	class Reflecting_1d
	{
	public:
		std::pair<double, double> boundaries;  // Lower and upper boundaries

    // Construct given lower and upper boundaries
		Reflecting_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

    // Check if position is out of bounds
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    // Enforce boundary condition
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += reflecting(state.position, boundaries);
      return true;
    }
	};

  // Reflecting boundaries along along each dimension
  // State must define: position (vector type)
	class Reflecting
	{
	public:
		const std::vector<std::pair<double, double>>
      boundaries;  // Lower and upper boundaries along each dimension

    // Construct given lower and upper boundaries along each dimension
		Reflecting(std::vector<std::pair<double, double>> boundaries)
		: boundaries(boundaries)
		{}

    // Check if position is out of bounds
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    // Enforce boundary condition
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
			for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
				state.position[dd] += reflecting(state.position[dd], boundaries[dd]);
      return true;
		}
	};

  // Reflecting boundary along dimension dd
  // State must define: position (vector type)
  template <std::size_t dd>
  class Reflecting_dim
  {
  public:
    
    std::pair<double, double> boundaries;  // Lower and upper boundaries along dimension dd
    static const std::size_t dim = dd;     // Dimension for outside visibility

    // Construct with boundaries at same distance from origin given domain halfwidth
    Reflecting_dim(double half_width)
    : boundaries{ -half_width, half_width }
    {}

    // Construct given lower and upper boundaries along dimension dim
    Reflecting_dim(std::pair<double, double> boundaries)
    : boundaries{ boundaries }
    {}

    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position[dd], boundaries); }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      state.position[1] +=
        boundary::reflecting(state.position[1], boundaries);
      return true;
    }
  };

  // Reflecting boundaries on the inside of an infinite cylinder in 3d
  // Cylinder longitudinal axis is dd_open, which may be 0 or 2
  // State must define: position (vector type)
  template <std::size_t dd_open = 0>
  class Open_RadialReflecting_3d
  {
    static_assert(dd_open != 0 && dd_open != 2,
                  "Cylinder axis dimension must be 0 or 2");
  public:
    const double radius;

    // Construct given domain radius along reflecting dimensions
    Open_RadialReflecting_3d(double radius)
    : radius{ radius }
    {}

    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      double radial_pos_sq = std::inner_product(
        position.cbegin()+begin_transverse, position.cbegin()+begin_transverse+1,
        position.cbegin()+begin_transverse, 0.);
      return radial_pos_sq > radius_sq;
    }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      if (!outOfBounds(state.position))
        return false;
      
      std::vector<double> position_radial_old(2);
      std::vector<double> jump_radial(2);
      std::vector<double> tangent_at_contact(2);
      
      for (std::size_t dd = 0; dd < position_radial_old.size(); ++dd)
      {
        jump_radial[dd] = state.position[dd+begin_transverse]
          -state_old.position[dd+begin_transverse];
        position_radial_old[dd] = state_old.position[dd+begin_transverse];
      }

      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump outside of boundary
        double pos_dot_jump = std::inner_product(
          position_radial_old.cbegin(), position_radial_old.cend(),
          jump_radial.cbegin(), 0.);
        double norm_sq_jump = std::inner_product(
          jump_radial.cbegin(), jump_radial.cend(),
          jump_radial.cbegin(), 0.);
        double radial_pos_sq_old = operation::abs_sq(position_radial_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux
          + (radius_sq-radial_pos_sq_old)/norm_sq_jump) - aux;

        // Place old position at contact
        for (std::size_t dd = 0; dd < position_radial_old.size(); ++dd)
          position_radial_old[dd] += alpha * jump_radial[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_radial_old[1]/radius;
        tangent_at_contact[1] = -position_radial_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = std::inner_product(
          tangent_at_contact.cbegin(), tangent_at_contact.cend(),
          jump_radial.cbegin(), 0.);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        double cos2 = std::cos(2.*contact_angle);
        double sin2 = std::sin(2.*contact_angle);

        // Outgoing jump at collision, rotate outer jump by -phi and flip the direction
        jump_radial[0] = (1-alpha)
          *(cos2*jump_radial[0] + sin2*jump_radial[1]);
        jump_radial[1] = (1-alpha)
          *(-sin2*jump_radial[0] + cos2*jump_radial[1]);
        
        // Place particle at reflected position
        for (std::size_t dd = 0; dd < position_radial_old.size(); ++dd)
          state.position[dd+begin_transverse] = position_radial_old[dd];
        
        // If inside boundary, done
        if(std::inner_product(state.position.cbegin()+begin_transverse,
                              state.position.cbegin()+begin_transverse+1,
                              state.position.cbegin()+begin_transverse, 0.))
          break;
      }
        
      return true;
    }

  private:
    const double radius_sq{ radius*radius };
    const std::size_t begin_transverse{ dd_open == 0 ? 1 : 0 };
  };
  
  // Place at closest void center if out of bounds
  // State must define: position (vector type)
  template <typename Grid, std::size_t dim>
  class ClosestVoidCenter
  {
  public:
    using KDTree =
      grid::KDTree_Grid_Mask<Grid, dim>;  // KDTree type to search for voids and solids
    
    // Construct given underlying grid, KDTree for void grid elements,
    // and KDTree for solid grid elements
    ClosestVoidCenter
    (Grid const& grid,
     KDTree const& kdtree_void,
     KDTree const& kdtree_solid)
    : grid{ grid }
    , kdtree_void{ kdtree_void }
    , kdtree_solid{ kdtree_solid }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      auto closest_void = kdtree_void.nearest_neighbor(position);
      auto closest_solid = kdtree_solid.nearest_neighbor(position);
      
      return closest_solid.second < closest_void.second;
    }
    
    // Check if out of bounds and provide closest void and closest solid
    // and respective squared distances to position
    template <typename Position>
    bool outOfBounds
    (Position const& position,
     std::pair<std::size_t, double>& closest_void,
     std::pair<std::size_t, double>& closest_solid) const
    {
      closest_void = kdtree_void.nearest_neighbor(position);
      closest_solid = kdtree_solid.nearest_neighbor(position);
      
      return closest_solid.second < closest_void.second;
    }
    
    // Check if out of bounds and provide closest void
    // and respective squared distance to position
    template <typename Position>
    bool outOfBounds
    (Position const& position,
     std::pair<std::size_t, double>& closest_void) const
    {
      closest_void = kdtree_void.nearest_neighbor(position);
      auto closest_solid = kdtree_solid.nearest_neighbor(position);
      
      return closest_solid.second < closest_void.second;
    }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      std::pair<std::size_t, double> closest_void;
      if (!outOfBounds(state.position, closest_void))
        return false;
      
      state.position =
        grid.cell_center(kdtree_void.points(closest_void.first));
      return true;
    }
    
  private:
    Grid const& grid;            // Underlying grid
    KDTree const& kdtree_void;   // KDTree for void grid cells
    KDTree const& kdtree_solid;  // KDTree for solid grid cells
  };
  
  // Place at closest void center if out of bounds
  // Switch to tumble state if boundary enforced (state.state = 2)
  // State must define: position (vector type)
  //                    state (integer type)
  template <typename Grid, std::size_t dim>
  class ClosestVoidCenter_Tumble
  {
  public:
    using KDTree =
      typename ClosestVoidCenter<Grid, dim>::KDTree;  // KDTree type to search for voids and solids
    
    // Construct given underlying grid, KDTree for void grid elements,
    // and KDTree for solid grid elements
    ClosestVoidCenter_Tumble
    (Grid const& grid,
     KDTree const& kdtree_void,
     KDTree const& kdtree_solid)
    : boundary{ grid, kdtree_void, kdtree_solid }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return boundary.outOfBounds(position); }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      // Check boundaries and enforce if needed
      if(!boundary(state, state_old))
        return false;
      
      // Set the tumble state
      state.state = 2;
      return true;
    }
    
  private:
    ClosestVoidCenter<Grid, dim> boundary;  // Boundary to place at closest void center
  };
  
  // Reflect off beads in a beadpack
  // State must define: position (vector type)
  template <typename BeadPack>
  class ReflectingBeads
  {
  public:
    // Construct given beadpack
    ReflectingBeads
    (BeadPack const& bead_pack)
    : bead_pack{ bead_pack }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return bead_pack.inside(position).first; }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      // If current position is inside a bead, apply reflection on that bead
      // Repeat until within periodic bounds and outside beads
      
      auto& position_new = state.position;
      
      auto inside_bead = bead_pack.inside(position_new);
      if (!inside_bead.first)
        return false;
      
      // Place current position at reflected position
      // and old position at reflection point
      // to compute displacement for further reflections if needed
      auto position_old = state_old.position;
      while (inside_bead.first)
      {
        geometry::ReflectOffSphere_OutsideToInside(position_new, position_old,
                                                   bead_pack.bead(inside_bead.second),
                                                   position_new, position_old);
        
        inside_bead = bead_pack.inside(position_new);
      }
      
      return true;
    }
    
  private:
    BeadPack const& bead_pack;  // Beadpack with beads to reflect off
  };
  
  // Reflect off beads in a beadpack
  // and enforce given periodic boundary conditions
  // State must define: position (vector type)
  template <typename BeadPack, typename PeriodicBoundary>
  class ReflectingBeads_Periodic
  {
  public:
    // Construct given beadpack object and boundary object to enforce periodic boundary
    ReflectingBeads_Periodic
    (BeadPack const& bead_pack, PeriodicBoundary const& periodic_boundary)
    : bead_pack{ bead_pack }
    , periodic_boundary{ periodic_boundary }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return
      bead_pack.inside(position).first
      || periodic_boundary.outOfBounds(position);
    }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      // Apply periodic boundary condition
      // If current position is inside a bead, apply reflection on that bead
      // Repeat until within periodic bounds and outside beads
      
      auto position_old = state_old.position;
      
      bool outofbounds_reflection;
      bool outofbounds_periodic;
      int counter_boundary_effects = 0;
      do
      {
        outofbounds_reflection = 0;
        outofbounds_periodic = 0;
        
        // Compute last displacement
        auto jump = operation::minus(state.position, position_old);
        outofbounds_periodic = periodic_boundary(state);
        // Adjust displacement if periodicity is enforced
        if (outofbounds_periodic)
          operation::minus(state.position, jump, position_old);
        
        // Reflect if inside bead
        // Place current position at reflected position
        // and old position at reflection point
        // to compute displacement for further reflections if needed
        auto inside_bead = bead_pack.inside(state.position);
        outofbounds_reflection = inside_bead.first;
        if (outofbounds_reflection)
          geometry::ReflectOffSphere_OutsideToInside(state.position, position_old,
                                                     bead_pack.bead(inside_bead.second),
                                                     state.position, position_old);
        
        // Check if anything happened, then end if not
        counter_boundary_effects += outofbounds_periodic + outofbounds_reflection;
      }
      while (outofbounds_periodic || outofbounds_reflection);
      
      return counter_boundary_effects;
    }
    
  private:
    BeadPack const& bead_pack;                  // Beadpack with beads to reflect off
    PeriodicBoundary const& periodic_boundary;  // Boundary object for periodic boundaries
  };
}

#endif /* Boundary_h */
