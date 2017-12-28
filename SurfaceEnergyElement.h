// -*- C++ -*-
#ifndef ELEMENT_SURFACE_ENERGY_CHEMOMECHANICAL
#define ELEMENT_SURFACE_ENERGY_CHEMOMECHANICAL

#include "Definitions.h"
#include "Utilities.h"
//#include "../../elements/Triangle.h"

namespace Elements {
namespace SurfaceGammaElement {

class Properties {
public:
  double _thickness;
	double _cmax;
	double _gamma_FePO4;
	double _gamma_LiFePO4;
	double _barrierCoeff;
  Properties(const double thickness, const double cmax, const double gamma_FePO4, const double gamma_LiFePO4, const double barrierCoeff) :
    _thickness(thickness), _cmax(cmax), _gamma_FePO4(gamma_FePO4), _gamma_LiFePO4(gamma_LiFePO4), _barrierCoeff(barrierCoeff) {
  }
};

template<unsigned QP>
class LinearTwoNodeSurfaceEnergyElement {

public:

  static const unsigned int     NumberOfNodes = 2;
  static const unsigned int			SpatialDimension = 2;
	//static const unsigned int			MechanicalDegreesOfFreedom = 3;
	static const unsigned int			TotalDegreesOfFreedom = 4;
	static const unsigned int			QuadPoints = QP;
  static const VTKCellType      VtkCellType = VTK_LINE;

  typedef SurfaceGammaElement::Properties												Properties;
	typedef Matrix<double, SpatialDimension, 1>										Point;
	typedef Matrix<double, TotalDegreesOfFreedom, 1>							TotalVariableVector;
	//typedef Matrix<double, MechanicalDegreesOfFreedom, 1>					MechanicalVector;
	//typedef array<MechanicalVector, NumberOfNodes>								MechanicalForces;
	typedef Matrix<double, SpatialDimension-1, 1>									QPoint;
typedef array<TotalVariableVector, NumberOfNodes>								Forces;
	typedef Matrix<double,
								 TotalDegreesOfFreedom*NumberOfNodes,
							   TotalDegreesOfFreedom*NumberOfNodes>						StiffnessMatrix;
	typedef NodeWithId<Point>																			Node;
  typedef array<TotalVariableVector, NumberOfNodes>							PrimitiveVariables;

  LinearTwoNodeSurfaceEnergyElement(
                          const array<Node, NumberOfNodes> & nodes,
                          const Properties & properties,
													const QuadratureRule<SpatialDimension-1, QP> * quadratureRule) :
													_properties(properties),
													_quadratureRule(quadratureRule){
														
		for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
      _nodeIds[nodeIndex] = nodes[nodeIndex]._id;
    }
														
		const Point & x0 = nodes[0]._position;
		const Point & x1 = nodes[1]._position;
	  const Point l = x0-x1;
		_halfLength = l.norm()/2.;
		_length = l.norm();
			
  }
	
	Forces
	computeForces(const PrimitiveVariables & primitives,
								const double phi,
								const double time) const {
		
		ignoreUnusedVariable(phi);
		ignoreUnusedVariable(time);
		Forces packedForces;
		packedForces[0].fill(0.);
		packedForces[1].fill(0.);
		array<QPoint,QP> qPoints = _quadratureRule->_points;

		Matrix<double, NumberOfNodes, 1> c;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
		}

		double delta_gamma = _properties._gamma_LiFePO4-_properties._gamma_FePO4;
		
		packedForces[0](2) = -(3*c(0)*c(0) + 2*c(0)*c(1) + c(1)*c(1) - 2*(2*c(0)+c(1))*_properties._cmax+ 0.06*_properties._cmax*_properties._cmax)*_length*delta_gamma*_properties._thickness/(2*_properties._cmax * _properties._cmax * _properties._cmax);
		
		packedForces[1](2) = -(3*c(1)*c(1) + 2*c(1)*c(0) + c(0)*c(0) - 2*(2*c(1)+c(0))*_properties._cmax+ 0.06*_properties._cmax*_properties._cmax)*_length*delta_gamma*_properties._thickness/(2*_properties._cmax*_properties._cmax*_properties._cmax);
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			
			packedForces[0](2) += _properties._barrierCoeff*(1./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) - 1./(cBarAtQuadPoint*cBarAtQuadPoint))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness;
		
			packedForces[1](2) += _properties._barrierCoeff*(1./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) - 1./(cBarAtQuadPoint*cBarAtQuadPoint))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness;
			
		}
		
		return packedForces;

	}
	
	StiffnessMatrix
	computeStiffnessMatrix(const PrimitiveVariables & primitives,
												 const double phi,
												 const double time) const {
		
		ignoreUnusedVariable(time);
		ignoreUnusedVariable(phi);
		StiffnessMatrix stiffness;
		stiffness.fill(0.);
		array<QPoint,QP> qPoints = _quadratureRule->_points;

		Matrix<double, NumberOfNodes, 1> c;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
		}
		
		double delta_gamma = _properties._gamma_LiFePO4-_properties._gamma_FePO4;
		
		stiffness(2,2) = -(6*c(0)+2*c(1)-4*_properties._cmax)*_length*_properties._thickness*delta_gamma/(2*_properties._cmax*_properties._cmax*_properties._cmax);
		
		stiffness(2,6) = -(2*c(0)+2*c(1)-2*_properties._cmax)*_length*_properties._thickness*delta_gamma/(2*_properties._cmax*_properties._cmax*_properties._cmax);
		
		stiffness(6,2) = -(2*c(1)+2*c(0)-2*_properties._cmax)*_length*_properties._thickness*delta_gamma/(2*_properties._cmax*_properties._cmax*_properties._cmax);
		
		stiffness(6,6) = -(6*c(1)+2*c(0)-4*_properties._cmax)*_length*_properties._thickness*delta_gamma/(2*_properties._cmax*_properties._cmax*_properties._cmax);
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			
			stiffness(2,2) += _properties._barrierCoeff*(2./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) + 2./(cBarAtQuadPoint*cBarAtQuadPoint*cBarAtQuadPoint))*((1.-lambda)/2.)*((1.-lambda)/2.)*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/_properties._cmax;
			
			stiffness(2,6) += _properties._barrierCoeff*(2./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) + 2./(cBarAtQuadPoint*cBarAtQuadPoint*cBarAtQuadPoint))*((1.-lambda)/2.)*((1.+lambda)/2.)*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/_properties._cmax;
			
			stiffness(6,2) +=_properties._barrierCoeff*(2./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) + 2./(cBarAtQuadPoint*cBarAtQuadPoint*cBarAtQuadPoint))*((1.+lambda)/2.)*((1.-lambda)/2.)*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/_properties._cmax;
			
			stiffness(6,6) += _properties._barrierCoeff*(2./((1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)*(1-cBarAtQuadPoint)) + 2./(cBarAtQuadPoint*cBarAtQuadPoint*cBarAtQuadPoint))*((1.+lambda)/2.)*((1.+lambda)/2.)*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/_properties._cmax;
			
		}
		
		return stiffness;
		
	}
	
  array<size_t, NumberOfNodes>
  getNodeIds() const {
    return _nodeIds;
  }

private:

  array<size_t, NumberOfNodes>                 _nodeIds;
  Properties                                   _properties;
	double																			 _halfLength;
	double																			 _length;
	const QuadratureRule<SpatialDimension-1, QP> * _quadratureRule;

};

}
}

#endif
