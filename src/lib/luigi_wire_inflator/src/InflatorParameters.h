#ifndef CELLPARAMETERS_H
#define CELLPARAMETERS_H

#include <utility>
#include <map>
#include <vector>

#include <vcg/space/point2.h>

struct TessellationParameters
{
	double max_area;
	double min_angle;

	TessellationParameters(void)
	{
		max_area = 0.0002;
		min_angle = 30;
	}
};

struct ParameterOperation
{
	typedef int NodeID;

	typedef enum
	{
		Translation,
		Radius
	} OperationType;

	OperationType                  type;

	// if type == Translation
	// nodes_displ maps the nodes to be translated to the respective translation vectors
	std::map<NodeID, vcg::Point2d> nodes_displ;

	// if type == Radius
	// nodes contains the nodes for which this radius-change operation is applied
	std::vector<NodeID>            nodes;
};

class CellParameters
{
public:
	CellParameters(void) : CellParameters(0) {;}

	CellParameters(size_t parameters_count)
	{
		m_parameters.resize(parameters_count);
	}

	virtual ~CellParameters(void)
	{
		;
	}

	/*!
	 * \brief numberOfParameters
	 * \return the number of pattern parameters
	 */
	virtual size_t numberOfParameters(void) const
	{
		return m_parameters.size();
	}

	/*!
	 * \brief parameter
	 * \param i the parameter index.
	 * \return the i-th parameter
	 */
	virtual double & parameter(int i)
	{
		assert(i>=0 && i<int(m_parameters.size()));
		return m_parameters[i];
	}

	virtual const double & cParameter(int i) const
	{
		assert(i>=0 && i<int(m_parameters.size()));
		return m_parameters[i];
	}

protected:
	std::vector<double> m_parameters;
};

#endif // CELLPARAMETERS_H
