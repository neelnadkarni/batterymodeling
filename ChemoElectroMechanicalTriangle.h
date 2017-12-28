// -*- C++ -*-
#ifndef ELEMENT_TRIANGLE_CHEMOMECHANICAL
#define ELEMENT_TRIANGLE_CHEMOMECHANICAL

#include "Definitions.h"
#include "Utilities.h"
//#include "../../elements/Triangle.h"

namespace Elements {
namespace TriangleForBatterySimulations {

class Properties {
public:
  double _thickness;
	double _cmax;
	
  Properties(const double thickness, const double cmax) :
    _thickness(thickness), _cmax(cmax) {
  }
};

	
template <class MaterialModel, unsigned QP>
class LinearChemoMechanical {

public:

  static const unsigned int     NumberOfNodes = 3;
  static const unsigned int			SpatialDimension = 2;
  // static const unsigned int			MechanicalDegreesOfFreedom = 3;
	// static const unsigned int			ChemicalDegreesOfFreedom = 1;
	static const unsigned int			PotentialDegreesOfFreedom = 1;
	static const unsigned int			TotalDegreesOfFreedom = 4;
	static const unsigned int     QuadPoints = QP;
  static const VTKCellType      VtkCellType = VTK_TRIANGLE;

  typedef TriangleForBatterySimulations::Properties							Properties;
	
  typedef Matrix<double, SpatialDimension, 1>										Point;
  // typedef Matrix<double, MechanicalDegreesOfFreedom, 1>					MechanicalVector;
	// typedef Matrix<double, ChemicalDegreesOfFreedom, 1>						ChemicalVector;
	typedef Matrix<double, PotentialDegreesOfFreedom, 1>					PotentialVector;
	typedef Matrix<double, TotalDegreesOfFreedom, 1>							TotalVariableVector;
	typedef array<TotalVariableVector, NumberOfNodes>							Forces;
	// typedef array<ChemicalVector, NumberOfNodes>									ChemicalForces;
	// typedef array<PotentialVector, NumberOfNodes>									PotentialForces;
	// typedef array<PotentialVector, NumberOfNodes>									ChemicalFluxForces;
  typedef NodeWithId<Point>																			Node;
  typedef array<TotalVariableVector, NumberOfNodes>							PrimitiveVariables;
  // typedef Matrix<double,
  //               NumberOfNodes*MechanicalDegreesOfFreedom,
  //               NumberOfNodes*MechanicalDegreesOfFreedom>			MechanicalStiffnessMatrix;
	typedef Matrix<double,
								 NumberOfNodes*TotalDegreesOfFreedom,
								 NumberOfNodes*TotalDegreesOfFreedom>						StiffnessMatrix;
	typedef Matrix<double,
                 NumberOfNodes*PotentialDegreesOfFreedom,
                 NumberOfNodes*PotentialDegreesOfFreedom>				NTNMatrix;
  typedef typename MaterialModel::Strain												Strain;
  typedef typename MaterialModel::Stress												Stress;
  typedef typename MaterialModel::TangentMatrix									TangentMatrix;

  double
  computeShapeFunctionDerivative(const unsigned int variable,
                                 const unsigned int index, const Point point) {
    ignoreUnusedVariables(point);
    switch (variable) {
    case 0:
      switch (index) {
      case 0:
        return 1;
      case 1:
        return 0;
      case 2:
        return -1;
      }
      break;
    case 1:
      switch (index) {
      case 0:
        return 0;
      case 1:
        return 1;
      case 2:
        return -1;
      }
    }
    printf("Cannot evaluate variable %u, "
           "shape function index %u\n", variable, index);
    exit(1);
  }

  LinearChemoMechanical(
                          const array<Node, NumberOfNodes> & nodes,
                          const Properties & properties,
                          const QuadratureRule<SpatialDimension, QP> * quadratureRule,
                          const MaterialModel * materialModel,
                          const array<typename MaterialModel::InternalVariables, QP> internalVariables = array<typename MaterialModel::InternalVariables, QP>()) :
    _quadratureRule(quadratureRule), _materialModel(materialModel),
    _properties(properties), _internalVariables(internalVariables){

    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
      _nodeIds[nodeIndex] = nodes[nodeIndex]._id;
    }

    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			
			const Point & qp = _quadratureRule->_points[qpIndex];
			
			Matrix<double, SpatialDimension, NumberOfNodes> shapeFunctionDerivatives;
			shapeFunctionDerivatives.fill(0);
			
			for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
        shapeFunctionDerivatives(0, nodeIndex) = computeShapeFunctionDerivative(0, nodeIndex, qp);
        shapeFunctionDerivatives(1, nodeIndex) = computeShapeFunctionDerivative(1, nodeIndex, qp);
      }
      //cout << "shapeFunctionDerivative: " << endl;
      //cout << shapeFunctionDerivatives << endl;

      Matrix<double, SpatialDimension, SpatialDimension> jacobian;
      jacobian.fill(0);
      for (size_t nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
        //cout << "nodIndex = " << nodeIndex << endl;
        for (unsigned int i = 0; i < SpatialDimension; ++i) {
          for (unsigned int j = 0; j < SpatialDimension; ++j) {
            jacobian(i, j) +=
              computeShapeFunctionDerivative(i, nodeIndex, qp)*nodes[nodeIndex]._position(j);
            //cout << "i,j = " << i << "," << j << endl;
            //cout << "position(j) = " << nodes[nodeIndex]._position(j) << endl;
          }
        }
      }
			
      Matrix<double, SpatialDimension, SpatialDimension> gamma = jacobian.inverse();
      //cout << " det gamma = " << gamma.determinant() << endl;
      _derivativeMatrices[qpIndex] = gamma * shapeFunctionDerivatives;
      //cout << "_derivativeMatrix" << endl;
      //cout << _derivativeMatrices[0]<< endl;
      //cout << "gamma:" <<endl;
      //cout << gamma << endl;

      _weightedJacobian[qpIndex] =
        jacobian.determinant() * _quadratureRule->_weights[qpIndex] * _properties._thickness;
    }
			
    double volume = 0.;
    for (unsigned int qpIndex = 0; qpIndex < QP; qpIndex++) {
      volume += _weightedJacobian[qpIndex];
    }
    _nodalWeights[0] = volume/NumberOfNodes;
    _nodalWeights[1] = volume/NumberOfNodes;
    _nodalWeights[2] = volume/NumberOfNodes;
			
    _ntnMatrix.fill(0.0);
    for (size_t nodeIndex1=0; nodeIndex1<NumberOfNodes*PotentialDegreesOfFreedom; nodeIndex1++) {
      for (size_t nodeIndex2=0; nodeIndex2<NumberOfNodes*PotentialDegreesOfFreedom; nodeIndex2++) {
				if(nodeIndex1==nodeIndex2){
					_ntnMatrix(nodeIndex1,nodeIndex2) = volume/6.;
				}
				else{
					_ntnMatrix(nodeIndex1,nodeIndex2) = volume/12.;
				}
      }
    }
			
		//	cout << _ntnMatrix << endl;
			
  }
	
	Forces
	computeForces(const PrimitiveVariables & primitives,
								const double phi,
								const double timeStep) const {
		
		ignoreUnusedVariable(phi);
		Forces packedForces;
		Matrix<double, TotalDegreesOfFreedom*NumberOfNodes, 1> forces;
		forces.fill(0.);
		const array<Strain, QP> strains = computeStrainsAtGaussPoints(primitives);
		array<Point,QP> qPoints = _quadratureRule->_points;
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			const Strain & strain = strains[qpIndex];
			Stress stress = _materialModel->computeStress(strain, _internalVariables[qpIndex], timeStep);
			
			double r = qPoints[qpIndex](0);
			double s = qPoints[qpIndex](1);
			
			Matrix<double, 1, PotentialDegreesOfFreedom*NumberOfNodes> nMatrix;
			nMatrix.fill(0.0);
			nMatrix(0,0) = r;   nMatrix(0,1) = s;   nMatrix(0,2) = 1.-r-s;
			
			Matrix<double, 8, TotalDegreesOfFreedom*NumberOfNodes> bMatrix; // 8 stress components
			bMatrix.fill(0.);
			for (unsigned int i=0; i<NumberOfNodes; i++) {
				bMatrix(0,4*i) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(1,4*i+1) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(2,4*i) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(2,4*i+1) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(3,4*i+2) = nMatrix(0,i);
				bMatrix(4,4*i+2) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(5,4*i+2) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(6,4*i+3) = -_derivativeMatrices[qpIndex](0,i)*timeStep;
				bMatrix(7,4*i+3) = -_derivativeMatrices[qpIndex](1,i)*timeStep;
			}
			
			forces += _weightedJacobian[qpIndex] * bMatrix.transpose() * stress;
			
		}
		
		//cout << _weightedJacobian[0] << endl;
		
		Matrix<double, NumberOfNodes, 1> muT;
		muT.fill(0.);
		muT(0) = primitives[0](3);
		muT(1) = primitives[1](3);
		muT(2) = primitives[2](3);
		
		Matrix<double, NumberOfNodes, 1> cT;
		cT.fill(0.);
		cT(0) = primitives[0](2);
		cT(1) = primitives[1](2);
		cT(2) = primitives[2](2);
		
		Matrix<double, NumberOfNodes, 1> muForce = _ntnMatrix*muT;
		Matrix<double, NumberOfNodes, 1> cForce = _ntnMatrix*cT;
		
		// cout << "The timestep is: "<< timeStep << endl << endl;
		
		forces(2) -= muForce(0);
		forces(6) -= muForce(1);
		forces(10) -= muForce(2);
		
		forces(3) += cForce(0);
		forces(7) += cForce(1);
		forces(11) += cForce(2);
		
		for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
			for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; dofIndex++) {
				packedForces[nodeIndex](dofIndex) = forces(TotalDegreesOfFreedom*nodeIndex+dofIndex);
			}
		}
		
		return packedForces;
	}

	StiffnessMatrix
	computeStiffnessMatrix(const PrimitiveVariables & primitives,
												 const double phi,
												 const double timeStep) const {
		
		ignoreUnusedVariable(phi);
		Matrix<double, NumberOfNodes*TotalDegreesOfFreedom,
		NumberOfNodes*TotalDegreesOfFreedom> stiffness;
		stiffness.fill(0.);
		const array<Strain, QP> strains = computeStrainsAtGaussPoints(primitives);
		array<Point,QP> qPoints = _quadratureRule->_points;
		for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
			const Strain & strain = strains[qpIndex];
			const TangentMatrix tangent = _materialModel->computeTangentMatrix(strain, _internalVariables[qpIndex], timeStep);
			double r = qPoints[qpIndex](0);
			double s = qPoints[qpIndex](1);
			
			Matrix<double, 1, PotentialDegreesOfFreedom*NumberOfNodes> nMatrix;
			nMatrix.fill(0.0);
			nMatrix(0,0) = r;   nMatrix(0,1) = s;   nMatrix(0,2) = 1.-r-s;
			
			Matrix<double, 8, TotalDegreesOfFreedom*NumberOfNodes> bMatrix; // 8 stress components
			Matrix<double, 8, TotalDegreesOfFreedom*NumberOfNodes> bTMatrix;
			bMatrix.fill(0.);
			bTMatrix.fill(0.);
			for (unsigned int i=0; i<NumberOfNodes; i++) {
				bMatrix(0,4*i) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(1,4*i+1) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(2,4*i) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(2,4*i+1) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(3,4*i+2) = nMatrix(0,i);
				bMatrix(4,4*i+2) = _derivativeMatrices[qpIndex](0,i);
				bMatrix(5,4*i+2) = _derivativeMatrices[qpIndex](1,i);
				bMatrix(6,4*i+3) = -_derivativeMatrices[qpIndex](0,i)*timeStep;
				bMatrix(7,4*i+3) = -_derivativeMatrices[qpIndex](1,i)*timeStep;
				
				bTMatrix(0,4*i) = _derivativeMatrices[qpIndex](0,i);
				bTMatrix(1,4*i+1) = _derivativeMatrices[qpIndex](1,i);
				bTMatrix(2,4*i) = 0.5*_derivativeMatrices[qpIndex](1,i);
				bTMatrix(2,4*i+1) = 0.5*_derivativeMatrices[qpIndex](0,i);
				bTMatrix(3,4*i+2) = nMatrix(0,i);
				bTMatrix(4,4*i+2) = _derivativeMatrices[qpIndex](0,i);
				bTMatrix(5,4*i+2) = _derivativeMatrices[qpIndex](1,i);
				bTMatrix(6,4*i+3) = _derivativeMatrices[qpIndex](0,i);
				bTMatrix(7,4*i+3) = _derivativeMatrices[qpIndex](1,i);
				
			}
			
			stiffness += _weightedJacobian[qpIndex] * bMatrix.transpose() * tangent * bTMatrix;
		}
		
		for (unsigned int j=0; j<NumberOfNodes; j++) {
			stiffness(3,4*j+2) += _ntnMatrix(0,j);
			stiffness(7,4*j+2) += _ntnMatrix(1,j);
			stiffness(11,4*j+2) += _ntnMatrix(2,j);
			stiffness(2,4*j+3) -= _ntnMatrix(0,j);
			stiffness(6,4*j+3) -= _ntnMatrix(1,j);
			stiffness(10,4*j+3) -= _ntnMatrix(2,j);
		}
			 
    return stiffness;
	}
	
	Forces
	computeForcesForPreviousConcentration(const PrimitiveVariables & primitives,
																				const double phi,
																				const double timeStep) const {
		ignoreUnusedVariable(phi);
		Forces packedForces;
		Matrix<double, TotalDegreesOfFreedom*NumberOfNodes, 1> forces;
		forces.fill(0.);
		
		Matrix<double, NumberOfNodes, 1> cT;
		cT.fill(0.);
		cT(0) = primitives[0](2);
		cT(1) = primitives[1](2);
		cT(2) = primitives[2][2];
		
		Matrix<double, NumberOfNodes, 1> muT;
		muT.fill(0.);
		muT(0) = primitives[0](3);
		muT(1) = primitives[1](3);
		muT(2) = primitives[2](3);
		
		Matrix<double, NumberOfNodes, 1> cForce = _ntnMatrix*cT;
		Matrix<double, NumberOfNodes, 1> muForce = _ntnMatrix*muT;
		
		forces(3) += cForce(0);
		forces(7) += cForce(1);
		forces(11) += cForce(2);
		
		//forces(2) += 0.*muForce(0);
		//forces(6) += 0.*muForce(1);
		//forces(10) += 0.*muForce(2);
		
		for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
			for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; dofIndex++) {
				packedForces[nodeIndex](dofIndex) = forces(TotalDegreesOfFreedom*nodeIndex+dofIndex);
			}
		}
		
		return packedForces;
	}
	
	array<Stress, QP>
  computeStressesAtGaussPoints(const PrimitiveVariables & primitives,
                               const double time) const {
    array<Stress, QP> stresses;
    const array<Strain, QP> strains = computeStrainsAtGaussPoints(primitives);
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      //stresses[qpIndex] = strains[qpIndex];
       stresses[qpIndex] = 
         _materialModel->computeStress(strains[qpIndex],
                                       _internalVariables[qpIndex], time);  // change this to plot stresses or strains
    }
    return stresses;
  }

  NTNMatrix
  computeNTNMatrix() const {
    return _ntnMatrix;
  }
	
	double
	computeAverageConcentration(const PrimitiveVariables & primitives,
															const double time) const {
		ignoreUnusedVariables(time);
		double c_avg = 0.;
		for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
			c_avg += primitives[nodeIndex](2)/NumberOfNodes;
		}
		
		return c_avg;
	}

  array<size_t, NumberOfNodes>
  getNodeIds() const {
    return _nodeIds;
  }

  array<double, NumberOfNodes>
  computeNodalWeights() const {
    return _nodalWeights;
  }

  void
  updateInternalVariables(const PrimitiveVariables & primitives,
                          const double time) {
    ignoreUnusedVariables(primitives, time);
  }
	
  array<Strain, QP>
  computeStrainsAtGaussPoints(const PrimitiveVariables & primitives) const {

    Matrix<double, NumberOfNodes, 1> u1, u2, c, mu;
    for (unsigned int nodeIndex=0; nodeIndex<NumberOfNodes; nodeIndex++) {
      u1(nodeIndex) = primitives[nodeIndex](0);
      u2(nodeIndex) = primitives[nodeIndex](1);
      c(nodeIndex) = primitives[nodeIndex](2);
      mu(nodeIndex) = primitives[nodeIndex](3);
    }

    array<Strain, QP> strains;
    array<Point,QP> qPoints = _quadratureRule->_points;
    for (unsigned int qpIndex=0; qpIndex<QP; qpIndex++) {
      strains[qpIndex].fill(0.);

      // compute eps11, eps22, eps12
      Vector2d uDerivatives = _derivativeMatrices[qpIndex] * u1;
      Vector2d vDerivatives = _derivativeMatrices[qpIndex] * u2;
      strains[qpIndex](0) = uDerivatives(0);
      strains[qpIndex](1) = vDerivatives(1);
      strains[qpIndex](2) = (uDerivatives(1) + vDerivatives(0))/2.;

      // compute c
			double r = qPoints[qpIndex](0);
			double s = qPoints[qpIndex](1);
			strains[qpIndex](3) = r*c(0) + s*c(1) + (1.-r-s)*c(2);;

      // compute c,1, c,2
			Vector2d cDerivatives = _derivativeMatrices[qpIndex] * c;
			strains[qpIndex](4) = cDerivatives(0);
			strains[qpIndex](5) = cDerivatives(1);
			
			// compute mu,1, mu,2
			Vector2d muDerivatives = _derivativeMatrices[qpIndex] * mu;
			strains[qpIndex](6) = muDerivatives(0);
			strains[qpIndex](7) = muDerivatives(1);
			
		}

    return strains;
  }

private:

  array<size_t, NumberOfNodes>                 _nodeIds;
  const QuadratureRule<SpatialDimension, QP> * _quadratureRule;
  const MaterialModel *                        _materialModel;
  Properties                                   _properties;
  array<typename MaterialModel::InternalVariables, QP> _internalVariables;
  array<float, QP>                             _weightedJacobian;
  array<Matrix<double,
               SpatialDimension, NumberOfNodes>, QP>    _derivativeMatrices;
  array<double, NumberOfNodes>                 _nodalWeights;
  NTNMatrix																		 _ntnMatrix;

};

}
}

#endif
