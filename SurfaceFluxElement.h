// -*- C++ -*-
#ifndef ELEMENT_SURFACE_FLUX_CHEMOMECHANICAL
#define ELEMENT_SURFACE_FLUX_CHEMOMECHANICAL

#include "Definitions.h"
#include "Utilities.h"
//#include "MeshUtilities.h"
//#include "../../elements/Triangle.h"

namespace Elements {
namespace SurfaceFluxElement {

class Properties {
public:
  double _thickness;
	double _R; // gas constant
	double _T; // temperature
	double _F; // Faraday constant
	double _cmax; // max concentration
	double _k; // reaction rate
  Properties(const double thickness, const double R, const double T, const double F, const double cmax, const double k) :
    _thickness(thickness), _R(R), _T(T), _F(F), _cmax(cmax), _k(k) {
  }
};

template <unsigned QP>
class LinearTwoNodeSurfaceFluxElement {

public:

  static const unsigned int     NumberOfNodes = 2;
  static const unsigned int			SpatialDimension = 2;
  static const unsigned int			PotentialDegreesOfFreedom = 1;
	static const unsigned int			TotalDegreesOfFreedom = 4; // u1, u2, c, mu
	static const unsigned int     QuadPoints = QP;
  static const VTKCellType      VtkCellType = VTK_LINE;

  typedef SurfaceFluxElement::Properties												Properties;
  typedef Matrix<double, SpatialDimension, 1>										Point;
	typedef Matrix<double, SpatialDimension-1, 1>									QPoint;
	typedef Matrix<double, PotentialDegreesOfFreedom, 1>					PotentialVector;
	typedef Matrix<double, TotalDegreesOfFreedom, 1>							TotalVariableVector;
	//typedef array<PotentialVector, NumberOfNodes>									ChemicalFluxForces;
	typedef array<TotalVariableVector, NumberOfNodes>							Forces;
	typedef Matrix<double,
								 TotalDegreesOfFreedom*NumberOfNodes,
								 TotalDegreesOfFreedom*NumberOfNodes>						StiffnessMatrix;
  typedef NodeWithId<Point>																			Node;
  typedef array<TotalVariableVector, NumberOfNodes>							PrimitiveVariables;
	typedef Matrix<double, 2 ,1>																	CurrentContributionMatrix;
	
  LinearTwoNodeSurfaceFluxElement(
                          const array<Node, NumberOfNodes> & nodes,
                          const Properties & properties,
                          const QuadratureRule<SpatialDimension-1, QP> * quadratureRule /* surface element*/) :
    _quadratureRule(quadratureRule),
    _properties(properties){

    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
      _nodeIds[nodeIndex] = nodes[nodeIndex]._id;
    }

		const Point & x0 = nodes[0]._position;
		const Point & x1 = nodes[1]._position;
	  const Point & l = x0-x1;
		_halfLength = l.norm()/2.;
			
  }
	
	
	Forces
	computeForces(const PrimitiveVariables & primitives,
								const double phi, // applied potential
								const double timeStep) const {
		
		ignoreUnusedVariable(time);
		Forces forces;
		forces[0].fill(0.);
		forces[1].fill(0.);
		array<QPoint,QP> qPoints = _quadratureRule->_points;
		
		Matrix<double, NumberOfNodes, 1> c, mu;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
			mu(nodeIndex) = primitives[nodeIndex](3);
		}
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			double j0 = 3*_properties._k*(1-cBarAtQuadPoint)*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
			
			forces[0](3) += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			forces[1](3) += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;

		}
		
		return forces;
		
	}
	
	StiffnessMatrix
	computeStiffnessMatrix(const PrimitiveVariables & primitives,
												 const double phi,
												 const double timeStep) const {
		ignoreUnusedVariable(time);
		StiffnessMatrix stiffness;
		stiffness.fill(0.);
		
		array<QPoint,QP> qPoints = _quadratureRule->_points;
		
		Matrix<double, NumberOfNodes, 1> c, mu;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
			mu(nodeIndex) = primitives[nodeIndex](3);
		}
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			double j0 = 3*_properties._k*(1-cBarAtQuadPoint)*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
			double dj0dc1 = 3*_properties._k*(1./_properties._cmax)*(-sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)) + 1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint)*(1.-cBarAtQuadPoint))*(-lambda+1.)/2.;
			double dj0dc2 = 3*_properties._k*(1./_properties._cmax)*(-sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)) + 1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint)*(1.-cBarAtQuadPoint))*(lambda+1.)/2.;
			
			stiffness(3,2) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			stiffness(3,3) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
			
			stiffness(3,6) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			stiffness(3,7) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
			
			stiffness(7,2) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			stiffness(7,3) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
			
			stiffness(7,6) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			stiffness(7,7) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
		}
		
		return stiffness;
		
	}
	
	double
	computeForceForCurrentEquation(const PrimitiveVariables & primitives,
																 const double phi,
																 const double time) {
		
		ignoreUnusedVariable(time);
		double force = 0.;
		array<QPoint,QP> qPoints = _quadratureRule->_points;
		
		Matrix<double, NumberOfNodes, 1> c, mu;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
			mu(nodeIndex) = primitives[nodeIndex](3);
		}
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			double j0 = 3*_properties._k*(1-cBarAtQuadPoint)*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
			
			force += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness;
		}
		
		return force;
		
	}
	
	// gradient of phi gradient of total current
	// c1, c2, mu1, mu2, and phi for the constant current boundary condition
	Matrix<double, 5, 1>
	computeStiffnessMatrixComponentsForCurrentEquation(const PrimitiveVariables & primitives,
																											const double phi,
																											const double time) const {
		ignoreUnusedVariable(time);
		Matrix<double, 5, 1> stiffness;
		stiffness.fill(0.);
		array<QPoint,QP> qPoints = _quadratureRule->_points;
		
		Matrix<double, NumberOfNodes, 1> c, mu;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
			mu(nodeIndex) = primitives[nodeIndex](3);
		}
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			
			double j0 = 3*_properties._k*(1-cBarAtQuadPoint)*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
			double dj0dc1 = 3*_properties._k*(1./_properties._cmax)*(-sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)) + 1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint)*(1.-cBarAtQuadPoint))*(-lambda+1.)/2.;
			double dj0dc2 = 3*_properties._k*(1./_properties._cmax)*(-sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)) + 1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint)*(1.-cBarAtQuadPoint))*(lambda+1.)/2.;
			
			
			stiffness(0) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T); // dmu1
			stiffness(1) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T); // dmu2
			stiffness(2) += -(2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T))*_properties._F; // dphi
			stiffness(3) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness; // dc1
			stiffness(4) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness; // dc2
		}
	
		return stiffness;
		
	}
	
	// last column components
	Matrix<double, 2, 1>
	computePhiStiffness(const PrimitiveVariables & primitives,
											const double phi,
											const double timeStep) {
		
		Matrix<double, 2, 1> stiffness;
		stiffness.fill(0.);
		array<QPoint,QP> qPoints = _quadratureRule->_points;
		
		Matrix<double, NumberOfNodes, 1> c, mu;
		for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
			c(nodeIndex) = primitives[nodeIndex](2);
			mu(nodeIndex) = primitives[nodeIndex](3);
		}
		
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			double lambda = qPoints[qpIndex](0);
			double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
			double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
			double j0 = 3*_properties._k*(1-cBarAtQuadPoint)*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
			
			stiffness(0) += 2*(-_properties._F)/(_properties._R*_properties._T)*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
			stiffness(1) += 2*(-_properties._F)/(_properties._R*_properties._T)*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
		}
		
		return stiffness;
		
	}
	
	array<size_t, NumberOfNodes>
	getNodeIds() const {
		return _nodeIds;
	}
	
private:

  array<size_t, NumberOfNodes>                 _nodeIds;
  const QuadratureRule<SpatialDimension-1, QP> * _quadratureRule;
  Properties                                   _properties;
	double _halfLength;
	
};

}
	
namespace SurfaceFluxElementDoyle {
		
		class Properties {
		public:
			double _thickness;
			double _R; // gas constant
			double _T; // temperature
			double _F; // Faraday constant
			double _cmax; // max concentration
			double _k; // reaction rate
			Properties(const double thickness, const double R, const double T, const double F, const double cmax, const double k) :
			_thickness(thickness), _R(R), _T(T), _F(F), _cmax(cmax), _k(k) {
			}
		};
		
		template <unsigned QP>
		class LinearTwoNodeSurfaceFluxElement {
			
		public:
			
			static const unsigned int     NumberOfNodes = 2;
			static const unsigned int			SpatialDimension = 2;
			static const unsigned int			PotentialDegreesOfFreedom = 1;
			static const unsigned int			TotalDegreesOfFreedom = 4; // u1, u2, c, mu
			static const unsigned int     QuadPoints = QP;
			static const VTKCellType      VtkCellType = VTK_LINE;
			
			typedef SurfaceFluxElementDoyle::Properties										Properties;
			typedef Matrix<double, SpatialDimension, 1>										Point;
			typedef Matrix<double, SpatialDimension-1, 1>									QPoint;
			typedef Matrix<double, PotentialDegreesOfFreedom, 1>					PotentialVector;
			typedef Matrix<double, TotalDegreesOfFreedom, 1>							TotalVariableVector;
			//typedef array<PotentialVector, NumberOfNodes>									ChemicalFluxForces;
			typedef array<TotalVariableVector, NumberOfNodes>							Forces;
			typedef Matrix<double,
			TotalDegreesOfFreedom*NumberOfNodes,
			TotalDegreesOfFreedom*NumberOfNodes>						StiffnessMatrix;
			typedef NodeWithId<Point>																			Node;
			typedef array<TotalVariableVector, NumberOfNodes>							PrimitiveVariables;
			typedef Matrix<double, 2 ,1>																	CurrentContributionMatrix;
			
			LinearTwoNodeSurfaceFluxElement(
																			const array<Node, NumberOfNodes> & nodes,
																			const Properties & properties,
																			const QuadratureRule<SpatialDimension-1, QP> * quadratureRule /* surface element*/) :
			_quadratureRule(quadratureRule),
			_properties(properties){
				
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
			_nodeIds[nodeIndex] = nodes[nodeIndex]._id;
		}
				
				const Point & x0 = nodes[0]._position;
				const Point & x1 = nodes[1]._position;
				const Point & l = x0-x1;
				_halfLength = l.norm()/2.;
				
			}
			
			
			Forces
			computeForces(const PrimitiveVariables & primitives,
										const double phi, // applied potential
										const double timeStep) const {
				
				ignoreUnusedVariable(time);
				Forces forces;
				forces[0].fill(0.);
				forces[1].fill(0.);
				array<QPoint,QP> qPoints = _quadratureRule->_points;
				
				Matrix<double, NumberOfNodes, 1> c, mu;
				for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
					c(nodeIndex) = primitives[nodeIndex](2);
					mu(nodeIndex) = primitives[nodeIndex](3);
				}
				
				for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
					
					double lambda = qPoints[qpIndex](0);
					double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
					double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
					double j0 = _properties._k*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
					
					forces[0](3) += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					forces[1](3) += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					
				}
				
				return forces;
				
			}
			
			StiffnessMatrix
			computeStiffnessMatrix(const PrimitiveVariables & primitives,
														 const double phi,
														 const double timeStep) const {
				ignoreUnusedVariable(time);
				StiffnessMatrix stiffness;
				stiffness.fill(0.);
				
				array<QPoint,QP> qPoints = _quadratureRule->_points;
				
				Matrix<double, NumberOfNodes, 1> c, mu;
				for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
					c(nodeIndex) = primitives[nodeIndex](2);
					mu(nodeIndex) = primitives[nodeIndex](3);
				}
				
				for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
					
					double lambda = qPoints[qpIndex](0);
					double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
					double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
					double j0 = _properties._k*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
					double dj0dc1 = _properties._k*(1./_properties._cmax)*(1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint))*(-lambda+1.)/2.;
					double dj0dc2 = _properties._k*(1./_properties._cmax)*(1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint))*(lambda+1.)/2.;
					
					stiffness(3,2) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					stiffness(3,3) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
					
					stiffness(3,6) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					stiffness(3,7) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
					
					stiffness(7,2) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					stiffness(7,3) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
					
					stiffness(7,6) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					stiffness(7,7) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T)*timeStep;
				}
				
				return stiffness;
				
			}
			
			double
			computeForceForCurrentEquation(const PrimitiveVariables & primitives,
																		 const double phi,
																		 const double time) {
				
				ignoreUnusedVariable(time);
				double force = 0.;
				array<QPoint,QP> qPoints = _quadratureRule->_points;
				
				Matrix<double, NumberOfNodes, 1> c, mu;
				for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
					c(nodeIndex) = primitives[nodeIndex](2);
					mu(nodeIndex) = primitives[nodeIndex](3);
				}
				
				for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
					
					double lambda = qPoints[qpIndex](0);
					double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
					double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
					double j0 = _properties._k*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
					
					force += 2*j0*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness;
				}
				
				return force;
				
			}
			
			// gradient of phi gradient of total current
			// c1, c2, mu1, mu2, and phi for the constant current boundary condition
			Matrix<double, 5, 1>
			computeStiffnessMatrixComponentsForCurrentEquation(const PrimitiveVariables & primitives,
																												 const double phi,
																												 const double time) const {
				ignoreUnusedVariable(time);
				Matrix<double, 5, 1> stiffness;
				stiffness.fill(0.);
				array<QPoint,QP> qPoints = _quadratureRule->_points;
				
				Matrix<double, NumberOfNodes, 1> c, mu;
				for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
					c(nodeIndex) = primitives[nodeIndex](2);
					mu(nodeIndex) = primitives[nodeIndex](3);
				}
				
				for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
					
					double lambda = qPoints[qpIndex](0);
					double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
					double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
					
					double j0 = _properties._k*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
					double dj0dc1 = _properties._k*(1./_properties._cmax)*(1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint))*(-lambda+1.)/2.;
					double dj0dc2 = _properties._k*(1./_properties._cmax)*(1./(2.*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint)))*(1-2.*cBarAtQuadPoint))*(lambda+1.)/2.;
					
					
					stiffness(0) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T); // dmu1
					stiffness(1) += 2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T); // dmu2
					stiffness(2) += -(2*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness/(_properties._R*_properties._T))*_properties._F; // dphi
					stiffness(3) += 2*dj0dc1*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness; // dc1
					stiffness(4) += 2*dj0dc2*sinh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness; // dc2
				}
				
				return stiffness;
				
			}
			
			// last column components
			Matrix<double, 2, 1>
			computePhiStiffness(const PrimitiveVariables & primitives,
													const double phi,
													const double timeStep) {
				
				Matrix<double, 2, 1> stiffness;
				stiffness.fill(0.);
				array<QPoint,QP> qPoints = _quadratureRule->_points;
				
				Matrix<double, NumberOfNodes, 1> c, mu;
				for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
					c(nodeIndex) = primitives[nodeIndex](2);
					mu(nodeIndex) = primitives[nodeIndex](3);
				}
				
				for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
					
					double lambda = qPoints[qpIndex](0);
					double muAtQuadPoint = (-lambda+1.)*mu(0)/2. + (lambda+1.)*mu(1)/2.;
					double cBarAtQuadPoint = ((-lambda+1.)*c(0)/2. + (lambda+1.)*c(1)/2.)/_properties._cmax;
					double j0 = _properties._k*sqrt(cBarAtQuadPoint*(1-cBarAtQuadPoint));
					
					stiffness(0) += 2*(-_properties._F)/(_properties._R*_properties._T)*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.-lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
					stiffness(1) += 2*(-_properties._F)/(_properties._R*_properties._T)*j0*cosh((muAtQuadPoint-_properties._F*phi)/(_properties._R*_properties._T))*(1.+lambda)/2.*_quadratureRule->_weights[qpIndex]*_halfLength*_properties._thickness*timeStep;
				}
				
				return stiffness;
				
			}
			
			array<size_t, NumberOfNodes>
			getNodeIds() const {
				return _nodeIds;
			}
			
		private:
			
			array<size_t, NumberOfNodes>                 _nodeIds;
			const QuadratureRule<SpatialDimension-1, QP> * _quadratureRule;
			Properties                                   _properties;
			double _halfLength;
			
		};
		
}
	
}

#endif
