#ifndef PATTERN_H
#define PATTERN_H

#include <vcg/complex/complex.h>

#include "clipperHelper.h"

template <class VcgMesh>
class Pattern2D
{
public:
	typedef typename VcgMesh::ScalarType     ScalarType;
	typedef typename vcg::Point2<ScalarType> Coord2Type;

	virtual ~Pattern2D(void) {}

	virtual bool generate(void) { return false; }
	virtual bool tessellate( VcgMesh & mesh) = 0;

	virtual const ClipperLib::Paths & getPaths(void)
	{
		return m_paths;
	}

	static const int ScaleFactor = (1 << 24);

	static ClipperLib::IntPoint convertToIntPoint(const Coord2Type & p)
	{
		return scaleToIntPoint<ScalarType>(p, ScaleFactor);
	}

	static Coord2Type convertToCoord2(const ClipperLib::IntPoint & p)
	{
		return (Coord2Type(p.X, p.Y) / ScaleFactor);
	}

protected:
	ClipperLib::Paths m_base_paths;
	ClipperLib::Paths m_paths;
};

#endif // PATTERN_H
