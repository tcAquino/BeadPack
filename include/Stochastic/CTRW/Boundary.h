//
//  Boundary.h
//  EColi
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

namespace boundary
{
	double periodic(double position, std::pair<double, double> const& boundaries)
	{
		double box_size = boundaries.second - boundaries.first;
		double pos = position - boundaries.first;
		return -std::floor(pos/box_size)*box_size;
	}
  
  std::pair<double, int> periodic_with_outside_info
  (double position, std::pair<double, double> const& boundaries)
  {
    double box_size = boundaries.second - boundaries.first;
    double pos = position - boundaries.first;
    int outside = int(std::floor(pos/box_size));
    return { -outside*box_size, outside };
  }

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

  bool outOfBounds_box(double position, std::pair<double,double> const& boundaries)
  {
    return position < boundaries.first || position > boundaries.second;
  }

  template <typename Position, typename Boundaries>
  bool outOfBounds_box(Position const& position, Boundaries const& boundaries)
  {
    position[0];
    for (std::size_t dd=0; dd<position.size(); ++dd)
      if (outOfBounds_box(position[dd], boundaries[dd]))
        return true;

    return false;
  }

	class DoNothing
	{
	public:
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return false; }

		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    { return false; }
	};

	class Periodic_1d
	{
	public:
		std::pair<double, double> boundaries;

		Periodic_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return OutOfBounds_Box(position, boundaries); }

		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += periodic(state.position, boundaries);
      return true;
    }
	};

	class Periodic
	{
	public:
		const std::vector<std::pair<double, double>> boundaries;

		Periodic(std::vector<std::pair<double, double>> boundaries)
		: boundaries(boundaries)
		{}

		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

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
  
  class Periodic_WithOutsideInfo
  {
  public:
    const std::vector<std::pair<double, double>> boundaries;

    Periodic_WithOutsideInfo(std::vector<std::pair<double, double>> boundaries)
    : boundaries(boundaries)
    {}

    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position, boundaries); }

    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
      {
        auto change_outside = periodic_with_outside_info(state.position[dd], boundaries[dd]);
        state.position[dd] += change_outside.first;
        state.periodicity[dd] += change_outside.second;
      }
      return true;
    }
  };

	class Reflecting_1d
	{
	public:
		std::pair<double, double> boundaries;

		Reflecting_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += reflecting(state.position, boundaries);
      return true;
    }
	};

	class Reflecting
	{
	public:
		const std::vector<std::pair<double, double>> boundaries;

		Reflecting(std::vector<std::pair<double, double>> boundaries)
		: boundaries(boundaries)
		{}

		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

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

  class Open_Reflecting_2d
  {
  public:
    std::pair<double, double> boundaries{};

    Open_Reflecting_2d(double half_width)
    : boundaries{ -half_width, half_width }
    {}

    Open_Reflecting_2d(std::pair<double, double> boundaries)
    : boundaries{ boundaries }
    {}

    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position[1], boundaries); }

    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      state.position[1] += boundary::reflecting(state.position[1], boundaries);
      return true;
    }
  };

  class Open_RadialReflecting_3d
  {
  public:
    const double radius;

    Open_RadialReflecting_3d(double radius)
    : radius{ radius }
    {}

    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      double radial_pos_sq = std::inner_product(std::next(position.cbegin()), position.cend(),
                                                std::next(position.cbegin()), 0.);
      return radial_pos_sq > radius_sq;
    }

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
        jump_radial[dd] = state.position[dd+1]-state_old.position[dd+1];
        position_radial_old[dd] = state_old.position[dd+1];
      }

      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump outside of boundary
        double pos_dot_jump = std::inner_product(position_radial_old.cbegin(), position_radial_old.cend(),
                                                 jump_radial.cbegin(), 0.);
        double norm_sq_jump = std::inner_product(jump_radial.cbegin(), jump_radial.cend(),
                                                 jump_radial.cbegin(), 0.);
        double radial_pos_sq_old = operation::abs_sq(position_radial_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux + (radius_sq-radial_pos_sq_old)/norm_sq_jump) - aux;

        // Place old position at contact
        for (std::size_t dd = 0; dd < position_radial_old.size(); ++dd)
          position_radial_old[dd] += alpha * jump_radial[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_radial_old[1]/radius;
        tangent_at_contact[1] = -position_radial_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = std::inner_product(tangent_at_contact.cbegin(), tangent_at_contact.cend(),
                                                     jump_radial.cbegin(), 0.);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        double cos2 = std::cos(2.*contact_angle);
        double sin2 = std::sin(2.*contact_angle);

        // Outgoing jump at collision, rotate outer jump by -phi and flip the direction
        jump_radial[0] = (1-alpha)*(cos2*jump_radial[0] + sin2*jump_radial[1]);
        jump_radial[1] = (1-alpha)*(-sin2*jump_radial[0] + cos2*jump_radial[1]);
        
        // Place particle at reflected position
        for (std::size_t dd = 0; dd < position_radial_old.size(); ++dd)
          state.position[dd+1] = position_radial_old[dd];
        
        // If inside boundary, done
        if(std::inner_product(std::next(state.position.cbegin()), state.position.cend(),
                              std::next(state.position.cbegin()), 0.))
          break;
      }
        
      return true;
    }

  private:
    const double radius_sq{ radius*radius };
  };
  
  template <typename Grid, std::size_t dim>
  class ClosestVoidCenter
  {
  public:
    using KDTree = grid::KDTree_Grid_Mask<Grid, dim>;
    
    ClosestVoidCenter
    (Grid const& grid,
     KDTree const& kdtree_void,
     KDTree const& kdtree_solid)
    : grid{ grid }
    , kdtree_void{ kdtree_void }
    , kdtree_solid{ kdtree_solid }
    {}
    
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      auto closest_void = kdtree_void.nearest_neighbor(position);
      auto closest_solid = kdtree_solid.nearest_neighbor(position);
      
      return closest_solid.second < closest_void.second;
    }
    
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
    
    template <typename Position>
    bool outOfBounds
    (Position const& position,
     std::pair<std::size_t, double>& closest_void) const
    {
      closest_void = kdtree_void.nearest_neighbor(position);
      auto closest_solid = kdtree_solid.nearest_neighbor(position);
      
      return closest_solid.second < closest_void.second;
    }
    
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      std::pair<std::size_t, double> closest_void;
      if (!outOfBounds(state.position, closest_void))
        return false;
      
      state.position = grid.cell_center(kdtree_void.points(closest_void.first));
      return true;
    }
    
  private:
    Grid const& grid;
    KDTree const& kdtree_void;
    KDTree const& kdtree_solid;
  };
  
  template <typename Grid, std::size_t dim>
  class ClosestVoidCenter_Tumble
  {
  public:
    using KDTree = typename ClosestVoidCenter<Grid, dim>::KDTree;
    
    ClosestVoidCenter_Tumble
    (Grid const& grid,
     KDTree const& kdtree_void,
     KDTree const& kdtree_solid)
    : boundary{ grid, kdtree_void, kdtree_solid }
    {}
    
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return boundary.outOfBounds(position); }
    
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if(!boundary(state, state_old))
        return false;
      
      state.state = 2;
      return true;
    }
    
  private:
    ClosestVoidCenter<Grid, dim> boundary;
  };
  
  template <typename BeadPack>
  class ReflectingBeads
  {
  public:
    ReflectingBeads
    (BeadPack const& bead_pack)
    : bead_pack{ bead_pack }
    {}
    
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return bead_pack.inside(position).first; }
    
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      auto& position_new = state.position;
      
      auto inside_bead = bead_pack.inside(position_new);
      if (!inside_bead.first)
        return false;
      
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
    BeadPack const& bead_pack;
  };
  
//  class PeriodicBccSymmetryPlanes
//  {
//  public:
//    const std::size_t dim = 3;
//
//    PeriodicBccSymmetryPlanes(double radius)
//    : radius{ radius }
//    , bcc_cell_side{ 2.*radius/std::sqrt(3.) }
//    , domain_size{ compute_domain_size() }
//    , normal_to_symmetry_planes{ compute_normal_to_symmetry_planes() }
//    , translation_symmetry_planes{ compute_translation_symmetry_planes() }
//    {}
//
//    template <typename Position>
//    bool outOfBounds(Position const& position)
//    {
//      std::vector<double> position_projected = project(position);
//
//      for (std::size_t dd = 0; dd < dim; ++dd)
//        if (position_projected[dd] > domain_size)
//          return true;
//
//      return false;
//    }
//
//    template <typename State>
//    bool operator()(State& state, State const& state_old = {})
//    {
//      auto& position = state.position;
//
//      std::vector<double> position_projected = project(position);
//
//      bool outside = 0;
//
//      std::vector<double> nr_domains_away;
//      for (std::size_t dd = 0; dd < dim; ++dd)
//      {
//        nr_domains_away.push_back(std::floor(position_projected[dd]/domain_size));
//        if (nr_domains_away.back() > 0.)
//          outside = 1;
//      }
//
//
//      for (std::size_t dd = 0; dd < dim; ++dd)
//        if (nr_domains_away[dd] > 0.)
//        {
//          auto translation = operation::times_scalar(nr_domains_away[dd], translation_symmetry_planes[dd]);
//          operation::minus_InPlace(position, translation);
//        }
//
//      return outside;
//    }
//
//  private:
//
//    template <typename Position>
//    auto project(Position const& position)
//    {
//      std::vector<double> position_projected;
//      for (std::size_t dd = 0; dd < dim; ++dd)
//        position_projected.push_back(operation::dot(position, normal_to_symmetry_planes[dd]));
//
//      return position_projected;
//    }
//
//    double compute_domain_size()
//    {
//      std::vector<double> x0 = { bcc_cell_side, 0., 0. };
//      std::vector<double> x1(3, bcc_cell_side);
//      std::vector<double> normal = { 0, 1., 1. };
//      operation::div_scalar_InPlace(normal, std::sqrt(2.));
//      return operation::dot(operation::minus(x1, x0), normal);
//    }
//
//    std::vector<std::vector<double>> compute_normal_to_symmetry_planes()
//    {
//      std::vector<std::vector<double>> normal_to_symmetry_planes{
//        { 0., 1., 1. },
//        { 1., 0., 1. },
//        { 1., 1., 0. }
//      };
//      for (auto& line : normal_to_symmetry_planes)
//        operation::div_scalar_InPlace(line, std::sqrt(2.));
//
//      return normal_to_symmetry_planes;
//    }
//
//    std::vector<std::vector<double>> compute_translation_symmetry_planes()
//    {
//      std::vector<std::vector<double>> translation_symmetry_planes = {
//        { -1., 1., 1. },
//        { 1., -1., 1. },
//        { 1., 1., -1. }
//      };
//      for (auto& line : translation_symmetry_planes)
//        operation::times_scalar_InPlace(bcc_cell_side, line);
//
//      return translation_symmetry_planes;
//    }
//
//    double radius;
//    double bcc_cell_side;
//    double domain_size;
//    std::vector<std::vector<double>> normal_to_symmetry_planes;
//    std::vector<std::vector<double>> translation_symmetry_planes;
//  };
  
  template <typename BeadPack, typename PeriodicBoundary>
  class ReflectingBeads_Periodic
  {
  public:
    ReflectingBeads_Periodic
    (BeadPack const& bead_pack, PeriodicBoundary const& periodic_boundary)
    : bead_pack{ bead_pack }
    , periodic_boundary{ periodic_boundary }
    {}
    
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return
      bead_pack.inside(position).first
      || periodic_boundary.outOfBounds(position);
    }
    
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      auto position_old = state_old.position;
      
      bool outofbounds_reflection;
      bool outofbounds_periodic;
      int counter_boundary_effects = 0;
      do
      {
        outofbounds_reflection = 0;
        outofbounds_periodic = 0;
        auto jump = operation::minus(state.position, position_old);
        outofbounds_periodic = periodic_boundary(state);
        if (outofbounds_periodic)
          operation::minus(state.position, jump, position_old);
        
        auto inside_bead = bead_pack.inside(state.position);
        outofbounds_reflection = inside_bead.first;
        if (outofbounds_reflection)
        {
          geometry::ReflectOffSphere_OutsideToInside(state.position, position_old,
                                                     bead_pack.bead(inside_bead.second),
                                                     state.position, position_old);
        }
        
        counter_boundary_effects += outofbounds_periodic + outofbounds_reflection;
      }
      while (outofbounds_periodic || outofbounds_reflection);
      
      return counter_boundary_effects;
    }
    
  private:
    BeadPack const& bead_pack;
    PeriodicBoundary const& periodic_boundary;
  };
}

#endif /* Boundary_h */
