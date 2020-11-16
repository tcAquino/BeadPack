//
//  MultiArray.h
//  CTRW_2
//
//	Adapted from https://stackoverflow.com/a/19725907/2684539
//  Created by Tomas Aquino on 9/28/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef MultiArray_h
#define MultiArray_h

#include <cassert>
#include <cstddef>
#include <vector>
#include <numeric>
#include <functional>
#include "general/Modular.h"

namespace useful
{
	template <typename value_t>
	class MultiArray
	{
	public:
    using value_type = value_t;
    
		explicit MultiArray(std::vector<std::size_t> dimensions)
		: dimensions(dimensions)
		, values(std::accumulate(dimensions.begin(), dimensions.end(), 1,
                             std::multiplies<std::size_t>()))
		{}

		explicit MultiArray(std::vector<std::size_t> dimensions, value_type value)
		: dimensions(dimensions)
		, values(std::accumulate(dimensions.begin(), dimensions.end(), 1,
                             std::multiplies<std::size_t>()), value)
		{}

		value_type const& operator[] (std::vector<std::size_t> const& indexes) const
		{
			return values[computeIndex(indexes)];
		}
		value_type& operator[] (std::vector<std::size_t> const& indexes)
		{
			return values[computeIndex(indexes)];
		}
		value_type const& operator[] (std::size_t index) const
		{
			return values[index];
		}
		value_type& operator[] (std::size_t index)
		{
			return values[index];
		}

		std::size_t computeIndex(std::vector<std::size_t> const& indexes) const
		{
      return modular::positionToIndex(indexes, dimensions);
		}

		std::vector<std::size_t> computeIndexes(std::size_t index) const
		{
      std::vector<std::size_t> indexes;
      computeIndexes(index, indexes);
      
			return indexes;
		}

		void computeIndexes(std::size_t index, std::vector<std::size_t>& indexes) const
		{
      modular::indexToPosition(index, dimensions, indexes);
		}

		std::size_t size() const
		{ return values.size(); }

		std::size_t size(std::size_t dd) const
		{
			return dimensions[ dd ];
		}

		std::size_t rank() const
		{ return dimensions.size(); }

		std::vector<std::size_t> const& sizes() const
		{ return dimensions; }

		auto begin()
		{ return values.begin(); }

		auto end()
		{ return values.end(); }

		auto begin() const
		{ return values.begin(); }

		auto end() const
		{ return values.end(); }

		auto cbegin() const
		{ return values.cbegin(); }

		auto cend() const
		{ return values.cend(); }

	private:
		std::vector<std::size_t> dimensions;
		std::vector<value_type> values;
	};
}


#endif /* MultiArray_h */
