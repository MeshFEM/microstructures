#ifndef CLIPPERHELPER_H
#define CLIPPERHELPER_H

#include <clipper.hpp>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <cmath>

namespace std
{
template <>
struct hash<ClipperLib::IntPoint>
{
	typedef ClipperLib::IntPoint argument_type;
	typedef size_t               result_type;

	result_type operator()(const argument_type& a) const
	{
		hash<ClipperLib::cInt> hasher;
		return hasher(a.X) * 31 + hasher(a.Y);
	}
};
}

/// IntPoint
inline const ClipperLib::cInt & clipperPointCCoord(const ClipperLib::IntPoint & p, int coord)
{
#ifdef use_xyz
	switch(coord)
	{
	case 0 : return p.X;
	case 1 : return p.Y;
	case 2 : return p.Z;
	}
#else
	switch(coord)
	{
	case 0 : return p.X;
	case 1 : return p.Y;
	}
#endif
	assert(0);
	return p.X;
}

inline ClipperLib::cInt & clipperPointCoord(ClipperLib::IntPoint & p, int coord)
{
#ifdef use_xyz
	switch(coord)
	{
	case 0 : return p.X;
	case 1 : return p.Y;
	case 2 : return p.Z;
	}
#else
	switch(coord)
	{
	case 0 : return p.X;
	case 1 : return p.Y;
	}
#endif
	assert(0);
	return p.X;
}

template <typename ScalarType>
inline ClipperLib::IntPoint scaleToIntPoint(const vcg::Point2<ScalarType> & a, ScalarType scaleF = 1)
{
	return ClipperLib::IntPoint(ClipperLib::cInt(double(a.X()) * scaleF), ClipperLib::cInt(double(a.Y()) * scaleF));
}

inline ClipperLib::IntPoint operator + (const ClipperLib::IntPoint & p, const ClipperLib::IntPoint & t)
{
	return ClipperLib::IntPoint(p.X + t.X, p.Y + t.Y);
}

inline ClipperLib::IntPoint operator - (const ClipperLib::IntPoint & p, const ClipperLib::IntPoint & t)
{
	return ClipperLib::IntPoint(p.X - t.X, p.Y - t.Y);
}

inline ClipperLib::IntPoint operator * (const ClipperLib::IntPoint & p, const double a)
{
	return ClipperLib::IntPoint(ClipperLib::cInt(p.X * a), ClipperLib::cInt(p.Y * a));
}

inline ClipperLib::IntPoint operator / (const ClipperLib::IntPoint & p, const double a)
{
	return ClipperLib::IntPoint(ClipperLib::cInt(p.X / a), ClipperLib::cInt(p.Y / a));
}

inline ClipperLib::IntPoint & operator += (ClipperLib::IntPoint & lhs, const ClipperLib::IntPoint & rhs)
{
	lhs.X += rhs.X;
	lhs.Y += rhs.Y;
	return lhs;
}

inline ClipperLib::IntPoint interpolate(const ClipperLib::IntPoint & a, const ClipperLib::IntPoint & b, const double t)
{
	vcg::Point2d p0(a.X, a.Y);
	vcg::Point2d p1(b.X, b.Y);
	vcg::Point2d p = (p0 * (1 - t)) + (p1 * t);
	ClipperLib::IntPoint ret = {ClipperLib::cInt(round(p.X())), ClipperLib::cInt(round(p.Y()))};
	return ret;
}

inline double distance(const ClipperLib::IntPoint & a, const ClipperLib::IntPoint & b)
{
	ClipperLib::IntPoint diff = (b-a);
	return vcg::math::Sqrt( double( (diff.X * diff.X) + (diff.Y * diff.Y) ) );
}

inline ClipperLib::IntPoint & operator *= (ClipperLib::IntPoint & lhs, const double rhs)
{
	lhs.X = ClipperLib::cInt(lhs.X * rhs);
	lhs.Y = ClipperLib::cInt(lhs.Y * rhs);
	return lhs;
}

inline bool operator < (const ClipperLib::IntPoint & lhs, const ClipperLib::IntPoint & rhs)
{
	ClipperLib::cInt xdiff = (lhs.X - rhs.X);
	if (xdiff < 0)
		return true;

	if (xdiff == 0)
		return (lhs.Y < rhs.Y);

	return false;
}

/// Path
inline ClipperLib::Path operator + (const ClipperLib::Path & p, const ClipperLib::IntPoint & t)
{
	ClipperLib::Path newP;
	newP.resize(p.size());
	for (size_t i=0; i<p.size(); ++i)
	{
		newP[i] = p[i] + t;
	}

	return newP;
}

inline ClipperLib::Path & operator *= (ClipperLib::Path & lhs, const double rhs)
{
	for (size_t i=0; i<lhs.size(); ++i)
		lhs[i] *= rhs;

	return lhs;
}

/// Paths
inline ClipperLib::Paths operator + (const ClipperLib::Paths & paths, const ClipperLib::IntPoint & t)
{
	ClipperLib::Paths ret;
	for (size_t i=0; i<paths.size(); ++i)
		ret.push_back((paths[i] + t));
	return ret;
}

inline ClipperLib::Paths & operator *= (ClipperLib::Paths & lhs, const double rhs)
{
	for (size_t i=0; i<lhs.size(); ++i)
		lhs[i] *= rhs;

	return lhs;
}

inline ClipperLib::Paths & operator << (ClipperLib::Paths & lhs, const ClipperLib::Paths & rhs)
{
	lhs.insert(lhs.end(), rhs.begin(), rhs.end());
	return lhs;
}

inline ClipperLib::IntRect getBounds(const ClipperLib::Paths & paths)
{
		ClipperLib::Clipper c;
		c.AddPaths(paths, ClipperLib::ptSubject, true);
		ClipperLib::IntRect rect = c.GetBounds();
		std::swap(rect.bottom, rect.top); // needed!! maybe due to clipper bug.
		return rect;
}

template <typename BoundsType>
inline bool isBoundary(const ClipperLib::IntPoint & p, const BoundsType & bounds);

template <>
inline bool isBoundary<ClipperLib::IntRect>(const ClipperLib::IntPoint & p, const ClipperLib::IntRect & bounds)
{
	return  (p.X == bounds.left   || p.X == bounds.right ||
	         p.Y == bounds.bottom || p.Y == bounds.top);
}

template <>
inline bool isBoundary<ClipperLib::Path>(const ClipperLib::IntPoint & p, const ClipperLib::Path & bounds)
{
	return (ClipperLib::PointInPolygon(p, bounds) <= 0);
}

inline ClipperLib::IntPoint getPointInHolePath(const ClipperLib::Path & poly)
{
	static const ClipperLib::cInt clipper_epsilon = 3;

	assert(poly.size() > 1);

	ClipperLib::IntPoint p0 = poly[0];
	ClipperLib::IntPoint p1 = poly[1];

	vcg::Point3d n(0, 0, 1);
	vcg::Point3d v(p1.X-p0.X, p1.Y-p0.Y, 0);
	v = v.normalized();

	vcg::Point3d outer_dir = v ^ n;

	vcg::Point2d midPoint(p0.X + (p1.X-p0.X) / 2, p0.Y + (p1.Y-p0.Y) / 2);
	midPoint += vcg::Point2d(outer_dir.X(), outer_dir.Y()) * clipper_epsilon;

	return ClipperLib::IntPoint(ClipperLib::cInt(round(midPoint.X())), ClipperLib::cInt(round(midPoint.Y())));
}

inline ClipperLib::Path getPointsInHoles(const ClipperLib::Paths & polygon)
{
	ClipperLib::Path ret;

	for (size_t i=0; i<polygon.size(); ++i)
	{
		const ClipperLib::Path & path = polygon[i];
		if (!ClipperLib::Orientation(path))
		{
			ret.push_back(getPointInHolePath(path));
		}
	}

	return ret;
}

inline ClipperLib::Path convertToPath(const ClipperLib::IntRect & rect)
{
	ClipperLib::Path path;
	path.push_back(ClipperLib::IntPoint(rect.left,  rect.bottom));
	path.push_back(ClipperLib::IntPoint(rect.right, rect.bottom));
	path.push_back(ClipperLib::IntPoint(rect.right, rect.top   ));
	path.push_back(ClipperLib::IntPoint(rect.left,  rect.top   ));
	return path;
}

#endif // CLIPPERHELPER_H
