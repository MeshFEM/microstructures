////////////////////////////////////////////////////////////////////////////////
// edge_hash.h
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Hash function for std::pair<int, int> (to be used for edges).
*/ 
//  Created:  05/24/2016 11:58:29
////////////////////////////////////////////////////////////////////////////////
#ifndef EDGE_HASH_H
#define EDGE_HASH_H

#include <functional>
#include <utility>

struct edge_hash {
	typedef std::pair<int, int> argument_type;
	typedef size_t result_type;

	result_type operator()(const argument_type& a) const
	{
		std::hash<int> hasher;
		argument_type e = a;
		if (e.first > e.second)
			std::swap(e.first, e.second);

		return hasher(e.first) * 31 + hasher(e.second);
	}
};

#endif /* end of include guard: EDGE_HASH_H */
