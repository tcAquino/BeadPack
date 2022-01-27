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
#include "general/Operations.h"
#include "Geometry/Shape.h"
#include "Geometry/SymmetryPlanes.h"
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
  template <typename Position = std::vector<double>,
  typename Boundaries = std::vector<std::pair<double, double>>>
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
		const std::pair<double, double> boundaries;  // Lower and upper boundaries

    // Construct given lower and upper boundaries
		Periodic_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

    // Check if position is out of bounds
		template <typename Position>
		bool outOfBounds(Position position) const
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
    
    template <typename Position, typename Projection = double>
    void translate(Position& position, Projection projection) const
    {
      operation::plus_InPlace(position,
                      (boundaries.second-boundaries.first)*projection);
    }
	};

  // Periodic boundaries along each dimension
  // State must define: position (vector type)
	class Periodic
	{
	public:
    // Lower and upper boundaries along each dimension
		const std::vector<std::pair<double, double>> boundaries;
    const std::vector<double> domain_dimensions;

    // Construct given lower and upper boundaries along each dimension
		Periodic(std::vector<std::pair<double, double>> boundaries)
    : boundaries{ boundaries }
    , domain_dimensions{ make_dimensions(boundaries) }
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
    
    // Translate position according to symmetry planes
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      operation::plus_InPlace(position,
                      operation::times(domain_dimensions, projections));
    }
    
  private:
    std::vector<double> make_dimensions
    (std::vector<std::pair<double, double>> const& boundaries) const
    {
      std::vector<double> dimensions;
      for (auto const& val : boundaries)
        dimensions.push_back(val.second - val.first);
        
      return dimensions;
    }
	};
  
  // Periodic boundary along dimension dd
  // State must define: position (vector type)
  template <std::size_t dd>
  class Periodic_dim
  {
  public:
    
    const std::pair<double, double> boundaries;  // Lower and upper boundaries along dimension dd
    static const std::size_t dim = dd;           // Dimension for outside visibility

    // Construct with boundaries at same distance from origin given domain halfwidth
    Periodic_dim(double half_width)
    : boundaries{ -half_width, half_width }
    {}

    // Construct given lower and upper boundaries along dimension dim
    Periodic_dim(std::pair<double, double> boundaries)
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
      state.position[dd] +=
        boundary::periodic(state.position[dd], boundaries);
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
    const std::vector<double> domain_dimensions;

    // Construct given lower and upper boundaries
    Periodic_WithOutsideInfo
    (std::vector<std::pair<double, double>> boundaries)
    : boundaries(boundaries)
    , domain_dimensions{ make_dimensions(boundaries) }
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
    
    // Translate position according to symmetry planes
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      operation::plus_InPlace(position,
                      operation::times(domain_dimensions, projections));
    }
    
  private:
    std::vector<double> make_dimensions
    (std::vector<std::pair<double, double>> const& boundaries) const
    {
      std::vector<double> dimensions;
      for (auto const& val : boundaries)
        dimensions.push_back(val.second - val.first);
        
      return dimensions;
    }
  };

  // Reflecting boundary in one dimension
  // State must define: position (scalar type)
	class Reflecting_1d
	{
	public:
		const std::pair<double, double> boundaries;  // Lower and upper boundaries

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
    
    const std::pair<double, double> boundaries;  // Lower and upper boundaries along dimension dd
    static const std::size_t dim = dd;           // Dimension for outside visibility

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
      state.position[dd] +=
        boundary::reflecting(state.position[dd], boundaries);
      return true;
    }
  };
  
  // Reflecting boundaries on the inside of a circle
  // State must define: position (vector type)
  class RadialReflecting_2d
  {
  public:
    const double radius;

    // Construct given domain radius along reflecting dimensions
    RadialReflecting_2d(double radius)
    : radius{ radius }
    {}

    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return operation::abs_sq(position) > radius_sq;
    }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      if (!outOfBounds(state.position))
        return false;
      
      auto position_old = state_old.position;
      std::vector<double> jump = operation::minus(state.position, position_old);
      std::vector<double> tangent_at_contact(2);
      
      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump up to the boundary
        double pos_dot_jump = operation::dot(position_old, jump);
        double norm_sq_jump = operation::abs_sq(jump);
        double pos_sq_old = operation::abs_sq(position_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux
                                 +(radius_sq-pos_sq_old)/norm_sq_jump)-aux;
        
        // Avoid numerical issues
        if (alpha <= 0. || useful::isnan(alpha))
        {
          operation::times_scalar_InPlace(radius/operation::abs(state.position),
                                          state.position);
          break;
        }
        
        // Place old position at contact
        for (std::size_t dd = 0; dd < 2; ++dd)
          position_old[dd] += alpha*jump[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_old[1]/radius;
        tangent_at_contact[1] = -position_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = operation::dot(tangent_at_contact, jump);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        // Outgoing jump at collision, rotate outer jump and flip the direction
        operation::rotate(jump, -2.*contact_angle);
        operation::times_scalar_InPlace(1.-alpha, jump);
        
        // Place particle at reflected position
        operation::plus(position_old, jump,
                        state.position);
        
        // If inside boundary, done
        if(!outOfBounds(state.position))
          break;
      }
      
      return true;
    }

  private:
    const double radius_sq{ radius*radius };
  };

  // Reflecting boundaries on the inside of an infinite cylinder in 3d
  // Cylinder longitudinal axis is dd_open, which may be 0 or 2
  // State must define: position (vector type)
  template <std::size_t dd_open = 0>
  class Open_RadialReflecting_3d
  {
    static_assert(dd_open == 0 || dd_open == 2,
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
        position.cbegin()+begin_transverse, position.cbegin()+begin_transverse+2,
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
      
      for (std::size_t dd = 0; dd < 2; ++dd)
      {
        jump_radial[dd] = state.position[dd+begin_transverse]
          -state_old.position[dd+begin_transverse];
        position_radial_old[dd] = state_old.position[dd+begin_transverse];
      }
      
      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump up to the boundary
        double pos_dot_jump = std::inner_product(
          position_radial_old.cbegin(), position_radial_old.cend(),
          jump_radial.cbegin(), 0.);
        double norm_sq_jump = std::inner_product(
          jump_radial.cbegin(), jump_radial.cend(),
          jump_radial.cbegin(), 0.);
        double radial_pos_sq_old = operation::abs_sq(position_radial_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux
                                 +(radius_sq-radial_pos_sq_old)/norm_sq_jump)-aux;
        
        // Avoid numerical issues
        if (alpha <= 0. || useful::isnan(alpha))
        {
           double radial_pos_sq = std::inner_product(
             state.position.cbegin()+begin_transverse,
             state.position.cbegin()+begin_transverse+2,
             state.position.cbegin()+begin_transverse, 0.);
           for (std::size_t dd = 0; dd < 2; ++dd)
             state.position[dd+begin_transverse] *= radius/std::sqrt(radial_pos_sq);
           break;
        }

        // Place old position at contact
        for (std::size_t dd = 0; dd < 2; ++dd)
          position_radial_old[dd] += alpha*jump_radial[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_radial_old[1]/radius;
        tangent_at_contact[1] = -position_radial_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = std::inner_product(
          tangent_at_contact.cbegin(), tangent_at_contact.cend(),
          jump_radial.cbegin(), 0.);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        // Outgoing jump at collision, rotate outer jump and flip the direction
        operation::rotate(jump_radial, -2.*contact_angle);
        operation::times_scalar_InPlace(1.-alpha, jump_radial);
        
        // Place particle at reflected position
        for (std::size_t dd = 0; dd < 2; ++dd)
          state.position[dd+begin_transverse] =
            position_radial_old[dd]+jump_radial[dd];
        
        // If inside boundary, done
        if(!outOfBounds(state.position))
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
        geometry::reflectOffSphere_outsideToInside(position_new, position_old,
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
          geometry::reflectOffSphere_outsideToInside(state.position, position_old,
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
  
  // Periodic boundaries along each symmetry plane
  // State must define: position (vector type)
  template <typename SymmetryPlanes>
  class Periodic_SymmetryPlanes
  {
  public:
    const SymmetryPlanes symmetry_planes;
    const double scale;
    const std::vector<double> origin;
    
    // Construct given symmetry plane object,
    // overall scale factor, and coordinate origin
    Periodic_SymmetryPlanes
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : symmetry_planes{ symmetry_planes }
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    // Construct given overall scale factor and coordinate origin
    Periodic_SymmetryPlanes
    (double scale = 1., std::vector<double> origin = {})
    : symmetry_planes{}
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}

    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
        if (std::floor(geometry::project(position, symmetry_planes, dd, scale, origin)) != 0.)
          return 1;
      return 0;
    }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      auto projections = place_in_unit_cell(state.position, symmetry_planes, scale, origin);
      
      for (auto const& val : projections)
      if (val)
        return true;
      
      return false;
    }
    
    // Translate position according to symmetry planes
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      geometry::translate(position, symmetry_planes, projections, scale);
    }
    
  };
  
  // Periodic boundaries along each symmetry plane
  // with information about where position would be outside domain
  // State must define: position (vector type)
  //                    periodicity (std::vector<int>) counting how many domains
  //                                                   have been traveled along each dimension
  template <typename SymmetryPlanes>
  class Periodic_SymmetryPlanes_WithOutsideInfo
  {
  public:
    const SymmetryPlanes symmetry_planes;
    const double scale;
    const std::vector<double> origin;

    // Construct given symmetry plane object,
    // overall scale factor, and coordinate origin
    Periodic_SymmetryPlanes_WithOutsideInfo
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : symmetry_planes{ symmetry_planes }
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    // Construct given overall scale factor and coordinate origin
    Periodic_SymmetryPlanes_WithOutsideInfo
    (double scale = 1., std::vector<double> origin = {})
    : symmetry_planes{}
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
        if (std::floor(geometry::project(position, symmetry_planes, dd, scale, origin)) != 0.)
          return 1;
      return 0;
    }

    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      auto projections =
        geometry::place_in_unit_cell(state.position, symmetry_planes, scale, origin);
      operation::plus_InPlace(state.periodicity, projections);
      
      for (auto const& val : projections)
        if (val != 0)
          return 1;
      
      return 0;
    }
    
    // Translate position according to symmetry planes
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      geometry::translate(position, symmetry_planes, projections, scale);
    }
  };
  
  // Reflect off beads in a beadpack with velocity info and acceleration field
  // State must define: position (vector type)
  //                    velocity (vector type)
  template <typename BeadPack, typename Acceleration>
  class ReflectingBeads_Velocity
  {
  public:
    // Construct given beadpack and background acceleration
    ReflectingBeads_Velocity
    (BeadPack const& bead_pack,
     Acceleration acceleration,
     double timestep,
     double radius = 0.,
     double restitution_coeff_norm = 1.,
     double restitution_coeff_tang = 1.)
    : bead_pack{ bead_pack }
    , acceleration{ acceleration }
    , timestep{ timestep }
    , radius{ radius }
    , restitution_coeff_norm{ restitution_coeff_norm }
    , restitution_coeff_tang{ restitution_coeff_tang }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return bead_pack.near(position, radius).first; }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      // If current position is inside a bead, apply reflection on that bead
      // Repeat until within periodic bounds and outside beads
      
      auto& position_new = state.position;
      auto& velocity_new = state.velocity;
      
      auto inside_bead = bead_pack.near(position_new);
      if (!inside_bead.first)
        return false;
      
      // Place current position at reflected position
      // and old position at reflection point
      // to compute displacement for further reflections if needed
      auto position_old = state_old.position;
      auto velocity_old = state_old.velocity;
      double time_to_contact = 0.;
      while (inside_bead.first)
      {
        geometry::reflectOffSphere_velocity_outsideToInside(
          position_new, position_old, velocity_old,
          bead_pack.bead(inside_bead.second),
          position_new, velocity_new,
          position_old, velocity_old, time_to_contact,
          timestep-time_to_contact,
          acceleration(state_old),
          restitution_coeff_norm, restitution_coeff_tang, radius);
        inside_bead = bead_pack.near(position_new, state_old.radius);
      }
      
      return true;
    }
    
  private:
    BeadPack const& bead_pack;  // Beadpack with beads to reflect off
    Acceleration acceleration; // Uniform background acceleration
    double timestep;
    double radius;
    double restitution_coeff_norm;
    double restitution_coeff_tang;
  };
 
  // Reflect off beads in a beadpack
  // and enforce reflecting boundary conditions on rectangular box
  // State must define: position (vector type)
  //                    velocity (vector type)
  template <typename BeadPack, typename Acceleration>
  class ReflectingBeads_Velocity_ReflectingBox
  {
  public:
    // Construct given beadpack object and boundary object to enforce reflecting boundary
    ReflectingBeads_Velocity_ReflectingBox
    (BeadPack const& bead_pack, Acceleration&& acceleration,
     std::vector<std::pair<double, double>> boundaries,
     double timestep,
     double radius = 0., double restitution_coeff_norm = 1.,
     double restitution_coeff_tang = 1.)
    : bead_pack{ bead_pack }
    , acceleration{ std::forward<Acceleration>(acceleration) }
    , timestep{ timestep }
    , radius{ radius }
    , restitution_coeff_norm{ restitution_coeff_norm }
    , restitution_coeff_tang{ restitution_coeff_tang }
    , boundaries{ boundaries }
    {
      for (auto& bound : this->boundaries)
      {
        bound.first += radius;
        bound.second -= radius;
      }
    }
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return
      bead_pack.near(position, radius).first
      || outOfBounds_box(position, boundaries);
    }
    
    // Enforce boundary condition
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      // Apply reflecting box boundary condition
      // If current position is inside a bead, apply reflection on that bead
      // Repeat until within periodic bounds and outside beads
      
      auto position_old = state_old.position;
      auto velocity_old = state_old.velocity;
      auto& position_new = state.position;
      auto& velocity_new = state.velocity;
      
      double time_to_contact = 0.;
      
      bool outofbounds_beadpack;
      bool outofbounds_reflecting;
      int counter_boundary_effects = 0;
      do
      {
        outofbounds_beadpack = 0;
        outofbounds_reflecting = 0;
        
        // Reflect if outside box
        // Place current position at reflected position
        // and old position at reflection point
        // to compute displacement for further reflections if needed
        outofbounds_reflecting = outOfBounds_box(position_new, boundaries);
        
        if (outofbounds_reflecting)
          reflect_box(position_new, position_old, velocity_old,
                      position_new, velocity_new,
                      position_old, velocity_old, time_to_contact,
                      acceleration(state_old));

        // Reflect if inside bead
        // Place current position at reflected position
        // and old position at reflection point
        // to compute displacement for further reflections if needed
        auto inside_bead = bead_pack.near(state.position, radius);
        outofbounds_beadpack = inside_bead.first;
        if (outofbounds_beadpack)
          geometry::reflectOffSphere_velocity_outsideToInside(
            position_new, position_old, velocity_old,
            bead_pack.bead(inside_bead.second),
            position_new, velocity_new,
            position_old, velocity_old, time_to_contact,
            timestep-time_to_contact,
            acceleration(state_old),
            restitution_coeff_norm, restitution_coeff_tang, radius);
        
        // Check if anything happened, then end if not
        counter_boundary_effects += outofbounds_beadpack + outofbounds_reflecting;
      }
      while (outofbounds_beadpack || outofbounds_reflecting);
      
      return counter_boundary_effects;
    }
    
    template <typename Position>
    void reflect_box
    (Position const& position_new, Position const& position_old,
     Position const& velocity_old,
     Position& position_reflected, Position& velocity_reflected,
     Position& position_contact, Position& velocity_contact,
     double& time_to_contact, Position const& acceleration) const
    {
      // Determine contact times from linear approximation of trajectory
      std::vector<double> contact_times;
      contact_times.reserve(position_new.size());
      // Compute time to contact along each dimension
      for (std::size_t dd = 0; dd < position_new.size(); ++dd)
      {
        if (position_new[dd] < boundaries[dd].first)
        {
          contact_times.push_back(std::min(
            timestep,
            (position_old[dd]-boundaries[dd].first)/velocity_old[dd]));
          continue;
        }
        if (position_new[dd] > boundaries[dd].second)
        {
          contact_times.push_back(
            (boundaries[dd].second-position_old[dd])/
            velocity_old[dd]);
        }
        contact_times.push_back(timestep);
      }
      auto first_contact = std::min_element(contact_times.cbegin(),
                                            contact_times.cend());
      time_to_contact = *first_contact;
      std::size_t dd_contact = first_contact-contact_times.begin();
      operation::linearOp(time_to_contact, velocity_old,
                          position_old,
                          position_contact);
      
      // Compute reflected velocity at contact
      operation::linearOp(time_to_contact, acceleration,
                          velocity_old,
                          velocity_contact);
      
      velocity_contact[dd_contact] = -velocity_contact[dd_contact];
      velocity_contact[dd_contact] *= restitution_coeff_norm;
      for (std::size_t dd = 0; dd < velocity_contact.size(); ++dd)
        if (dd != dd_contact)
          velocity_contact[dd] *= restitution_coeff_tang;
      
      // Compute reflected position and velocity
      double time_leftover = timestep-time_to_contact;
      operation::linearOp(time_leftover, velocity_contact,
                          position_contact,
                          position_reflected);
      operation::linearOp(time_leftover, acceleration,
                          velocity_contact,
                          velocity_reflected);
    }
    
  private:
    BeadPack const& bead_pack;
    Acceleration acceleration;
    double timestep;
    double radius;
    double restitution_coeff_norm;
    double restitution_coeff_tang;
    std::vector<std::pair<double, double>> boundaries;
  };
  template
  <typename BeadPack, typename Acceleration>
  ReflectingBeads_Velocity_ReflectingBox
  (BeadPack const&, Acceleration&&,
   std::vector<std::pair<double, double>>,
   double, double, double, double) ->
  ReflectingBeads_Velocity_ReflectingBox<BeadPack, Acceleration>;
}

#endif /* Boundary_h */
