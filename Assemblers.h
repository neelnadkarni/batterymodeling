// -*- C++ -*-
#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "Definitions.h"
#include "Utilities.h"
//#include "MeshUtilities.h"
#include <complex>
using namespace std;

namespace Assemblers {

namespace Utilities {

template <class Element>
void
updateInternalVariablesUtility(const vector<typename Element::TotalVariableVector> & primitives,
                               const double time,
                               vector<Element> * elements) {
  for (size_t elementIndex = 0; elementIndex < elements->size(); ++elementIndex) {
    Element & element = elements->at(elementIndex);
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::TotalVariableVector,
                Element::NumberOfNodes> elementPrimitives =
      ::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
                                                                  primitives);
    element.updateInternalVariables(elementPrimitives, time);

  }
}

	
template <class Element>
void
assembleForceVectorUtility(const vector<typename Element::TotalVariableVector> & allPrimitives,
													 const double phi,
													 const vector<Element> & elements,
													 const double time,
													 Eigen::VectorXd * forces) {
  /*const size_t numberOfDofs = allPrimitives.size() * Element::TotalDegreesOfFreedom+1; // 1 extra for phi
  if (forces->size() != numberOfDofs) {
		errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
									 "and numberOfDofs = %zu\n", forces->size(),
									 numberOfDofs);
		exit(1);
	}*/
		
  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		// get element displacements
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																																allPrimitives);
		
		try {
			// assemble local contribution
			typename Element::Forces elementForces =
			element.computeForces(elementPrimitives, phi, time);
			for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
					 ++nodeIndex) {
				size_t nodeId = elementNodeIds[nodeIndex];
				for (size_t i = 0; i < Element::TotalDegreesOfFreedom; ++i) {
					(*forces)(nodeId * Element::TotalDegreesOfFreedom + i) +=
					elementForces[nodeIndex](i);
				}
			}
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's force vector\n",
										 elementIndex);
			throw;
		}
	}
}
	
template <class Element>
void
assembleAverageConcentrationUtility(const vector<typename Element::TotalVariableVector>
																		allPrimitives,
																		const vector<Element> & elements,
																		const double time,
																		double * averageConcentration) {

  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		// get element displacements
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																														 allPrimitives);
		
		try {
			// assemble local contribution
			(*averageConcentration) +=
			element.computeAverageConcentration(elementPrimitives, time)/elements.size();
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's average concentration \n",
										 elementIndex);
			throw;
		}
	}
}

template <class Element>
void
assemblePreviousConcentrationForceVectorUtility(const vector<typename Element::TotalVariableVector> & allPrimitives,
																								const double phi,
																								const vector<Element> & elements,
																								const double timeStep,
																								Eigen::VectorXd * forces) {
  /*const size_t numberOfDofs = allPrimitives.size() * Element::TotalDegreesOfFreedom+1;
  if (forces->size() != numberOfDofs) {
		errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
									 "and numberOfDofs = %zu\n", forces->size(),
									 numberOfDofs);
		exit(1);
	}*/
		
  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		// get element displacements
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																														 allPrimitives);
		
		try {
			// assemble local contribution
			typename Element::Forces elementForces =
			element.computeForcesForPreviousConcentration(elementPrimitives, phi, timeStep);
			for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
					 ++nodeIndex) {
				size_t nodeId = elementNodeIds[nodeIndex];
				for (size_t i = 0; i < Element::TotalDegreesOfFreedom; ++i) {
					(*forces)(nodeId * Element::TotalDegreesOfFreedom + i) +=
					elementForces[nodeIndex](i);
				}
			}
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's force vector\n",
										 elementIndex);
			throw;
		}
	}
}
	
/*
template <class Element>
void
assembleChemicalForceVectorUtility(const vector<typename Element::TotalVariableVector> & allPrimitives,
																			 const vector<Element> & elements,
																			 const double time,
																			 Eigen::VectorXd * forces) {
  const size_t numberOfDofs = allPrimitives.size() * Element::ChemicalDegreesOfFreedom;
  if (forces->size() != numberOfDofs) {
		errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
									 "and numberOfDofs = %zu\n", forces->size(),
									 numberOfDofs);
		exit(1);
	}
		
  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		// get element displacements
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																														 allPrimitives);
		
		try {
			// assemble local contribution
			typename Element::ChemicalForces elementForces =
			element.computeChemicalForces(elementPrimitives, time);
			for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
					 ++nodeIndex) {
				size_t nodeId = elementNodeIds[nodeIndex];
				for (size_t i = 0; i < Element::ChemicalDegreesOfFreedom; ++i) {
					(*forces)(nodeId * Element::ChemicalDegreesOfFreedom + i) +=
					elementForces[nodeIndex](i);
				}
			}
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's force vector\n",
										 elementIndex);
			throw;
		}
	}
}*/

	
/*template <class Element>
void
assembleChemicalFluxForceVectorUtility(const vector<typename Element::TotalVariableVector> & allPrimitives,
																		 const vector<Element> & elements,
																		 const double time,
																		 const double phi,
																		 Eigen::VectorXd * forces) {
  const size_t numberOfDofs = allPrimitives.size() * Element::PotentialDegreesOfFreedom;
  if (forces->size() != numberOfDofs) {
		errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
									 "and numberOfDofs = %zu\n", forces->size(),
									 numberOfDofs);
		exit(1);
	}
		
  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		// get element displacements
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																														 allPrimitives);
		
		try {
			// assemble local contribution
			typename Element::ChemicalFluxForces elementForces =
			element.computeChemicalFluxForces(elementPrimitives, phi, time);
			for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
					 ++nodeIndex) {
				size_t nodeId = elementNodeIds[nodeIndex];
				for (size_t i = 0; i < Element::PotentialDegreesOfFreedom; ++i) {
					(*forces)(nodeId * Element::PotentialDegreesOfFreedom + i) +=
					elementForces[nodeIndex](i);
				}
			}
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's force vector\n",
										 elementIndex);
			throw;
		}
	}
}*/

template <class Element, class MatrixExtractor>
void
assembleMatrixUtility(const size_t numberOfNodes,
												const vector<Element> & elements,
												const MatrixExtractor & matrixExtractor,
												Eigen::SparseMatrix<double> * matrix) {
  /*const size_t numberOfDofs = numberOfNodes * Element::TotalDegreesOfFreedom+1; // 1 extra for phi
  // doesn't matter now for sparse matricies does it?
  if (matrix->rows() != numberOfDofs || matrix->cols() != numberOfDofs) {
		throwException("cannot assembleMatrix with matrix->rows() = %d "
									 "and matrix->rows() = %d and numberOfDofs = %zu,"
									 "numberOfNodes = %lu,"
									 "DegreesOfFreedom = %u\n",
									 matrix->rows(), matrix->cols(), numberOfDofs,
									 numberOfNodes,
									 Element::TotalDegreesOfFreedom);
	}*/
  std::vector<Eigen::Triplet<double> > tripletList;
  Eigen::SparseMatrix<double> tempStorageMatrix(matrix->rows(),matrix->cols());
  // esimate the size based on the number of non zero entries
  // in the first matrix
  if (elements.size() > 0){
		const Element & exampleElement = elements[0];
		const typename MatrixExtractor::MatrixType elementMatrixExample =
		matrixExtractor(exampleElement);
		// SO this could be better thought out but it's the best I have right now
		tripletList.reserve(elements.size() * 20);
	}
  for (size_t elementIndex = 0; elementIndex < elements.size();
			 ++elementIndex) {
		const Element & element = elements[elementIndex];
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		
		try {
			const typename MatrixExtractor::MatrixType elementMatrix =
			matrixExtractor(element);
			
			// assemble local contribution
			for (size_t nodeIndex1 = 0; nodeIndex1 < elementNodeIds.size();
					 ++nodeIndex1) {
				const size_t nodeId1 = elementNodeIds[nodeIndex1];
				for (size_t nodeIndex2 = 0; nodeIndex2 < elementNodeIds.size();
						 ++nodeIndex2) {
					const size_t nodeId2 = elementNodeIds[nodeIndex2];
					for (size_t i = 0; i < Element::TotalDegreesOfFreedom; ++i) {
						for (size_t j = 0; j < Element::TotalDegreesOfFreedom; ++j) {
							// heaven help me if this is wrong ... seriously
							tripletList.push_back(Eigen::Triplet<double>
																		(nodeId1 * Element::TotalDegreesOfFreedom + i,
																		 nodeId2 * Element::TotalDegreesOfFreedom + j,
																		 elementMatrix(nodeIndex1 *
																									 Element::TotalDegreesOfFreedom + i,
																									 nodeIndex2 *
																									 Element::TotalDegreesOfFreedom + j)));
						}
					}
				}
			}
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble element %zu's matrix\n",
										 elementIndex);
			throw;
		}
	}
  tempStorageMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
  (*matrix) += tempStorageMatrix;
}
	
	
	
template <class Element>
struct StiffnessMatrixExtractor {
	
  typedef typename Element::StiffnessMatrix MatrixType;
	
  const vector<typename Element::TotalVariableVector> & _primitives;
	const double _phi;
	const double _time;
	StiffnessMatrixExtractor(const vector<typename Element::TotalVariableVector> & primitives,
													 const double phi,
													 const double time) :
		_primitives(primitives), _phi(phi), _time(time) {
	}
		
  MatrixType
  operator()(const Element & element) const {
		const array<size_t, Element::NumberOfNodes> elementNodeIds =
		element.getNodeIds();
		const array<typename Element::TotalVariableVector, Element::NumberOfNodes>
		elementPrimitives =
		::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
																																_primitives);
		return element.computeStiffnessMatrix(elementPrimitives, _phi, _time);
	}
	private:
  StiffnessMatrixExtractor();
};
	
	
template <class Element>
void
assembleStiffnessMatrixUtility(const vector<typename Element::TotalVariableVector> & allPrimitives,
															 const double phi,
                               const vector<Element> & elements,
                               const double time,
                               Eigen::SparseMatrix<double> * matrix) {
  assembleMatrixUtility<Element,
                        StiffnessMatrixExtractor<Element> >(allPrimitives.size(),
                                                            elements,
                                                            StiffnessMatrixExtractor<Element>(allPrimitives,phi,time),
                                                            matrix);
}

template <class Element>
vector<typename Element::Stress>
computeElementStresses(const vector<typename Element::TotalVariableVector> & allPrimitives,
                       const vector<Element> & elements,
                       const double time) {

  vector<typename Element::Stress> allElementStresses(elements.size());
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<typename Element::TotalVariableVector, Element::NumberOfNodes> elementPrimitives =
      ::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
                                                                  allPrimitives);
    array<typename Element::Stress, Element::QuadPoints> elementStresses =
      element.computeStressesAtGaussPoints(elementPrimitives, time);
    typename Element::Stress average = Element::Stress::Zero();
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      average += elementStresses[qpIndex];
    }
    average /= Element::QuadPoints;
    allElementStresses[elementIndex] = average;
  }
  return allElementStresses;
}

template <class Element>
vector<typename Element::Stress>
computeNodalStresses(const vector<typename Element::TotalVariableVector> & allPrimitives,
                     const vector<Element> & elements,
                     const double time) {

  const size_t numberOfNodes = allPrimitives.size();

  vector<typename Element::Stress> nodalStresses(numberOfNodes,
                                                 Element::Stress::Zero());
  vector<double> volumeSums(numberOfNodes, 0);

  vector<typename Element::Stress> elementStresses =
    computeElementStresses(allPrimitives, elements, time);

  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<double, Element::NumberOfNodes> elementWeights =
      element.computeNodalWeights();
    const typename Element::Stress & elementStress = elementStresses[elementIndex];

    for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
      nodalStresses[elementNodeIds[nodeIndex]] +=
        elementStress / elementWeights[nodeIndex];
      volumeSums[elementNodeIds[nodeIndex]] += 1./elementWeights[nodeIndex];
    }
  }

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
    nodalStresses[nodeIndex] /= volumeSums[nodeIndex];
  }
  return nodalStresses;
}

template <class Element>
vector<typename Element::Strain>
computeElementStrains(const vector<typename Element::TotalVariableVector> & allPrimitives,
                      const vector<Element> & elements) {

  vector<typename Element::Strain> allElementStrains(elements.size());
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<typename Element::TotalVariableVector, Element::NumberOfNodes> elementPrimitives =
      ::Utilities::getElementPrimitivesFromGlobalList<Element>(elementNodeIds,
                                                                  allPrimitives);
    array<typename Element::Strain, Element::QuadPoints> elementStrains =
      element.computeStrainsAtGaussPoints(elementPrimitives);
    typename Element::Strain average = Element::Strain::Zero();
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      average += elementStrains[qpIndex];
    }
    average /= Element::QuadPoints;
    allElementStrains[elementIndex] = average;
  }
  return allElementStrains;
}

template <class Element>
vector<typename Element::Strain>
computeNodalStrains(const vector<typename Element::TotalVariableVector> & allPrimitives,
                    const vector<Element> & elements) {

  const size_t numberOfNodes = allPrimitives.size();

  vector<typename Element::Strain> nodalStrains(numberOfNodes,
                                                Element::Strain::Zero());
  vector<double> volumeSums(numberOfNodes, 0);

  vector<typename Element::Strain> elementStrains =
    computeElementStrains(allPrimitives, elements);

  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<double, Element::NumberOfNodes> elementWeights =
      element.computeNodalWeights();
    const typename Element::Strain & elementStrain = elementStrains[elementIndex];

    for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
      nodalStrains[elementNodeIds[nodeIndex]] +=
        elementStrain / elementWeights[nodeIndex];
      volumeSums[elementNodeIds[nodeIndex]] += 1./elementWeights[nodeIndex];
    }
  }

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    nodalStrains[nodeIndex] /= volumeSums[nodeIndex];
  }
  return nodalStrains;
}

}

/*
template <class PhysicalElementType,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType> >
class AssemblerWithOrderParameters : public Assembler<PhysicalElementType, ExternalForceAssembler> {

  typedef Assembler<PhysicalElementType, ExternalForceAssembler> Base;
public:
  typedef typename PhysicalElementType::Vector          ElementVector;
  typedef typename PhysicalElementType::OrderParameter  OrderParameter;
  static const unsigned int      OrderParameterDofs   = PhysicalElementType::OrderParameterDofs;

  AssemblerWithOrderParameters(const size_t numberOfNodes) : Base(numberOfNodes) {
  }

  template <class QuadratureRule, class MaterialModel>
  AssemblerWithOrderParameters(const SingleElementMesh<PhysicalElementType> & mesh,
                               const typename PhysicalElementType::Properties & elementProperties,
                               const QuadratureRule & quadratureRule,
                               const MaterialModel & materialModel) :
    Base(mesh, elementProperties, quadratureRule, materialModel) {
  }

  void
  updateOrderParameters(const vector<OrderParameter> & orderParameters)  {
    Assemblers::Utilities::updateOrderParametersUtility(orderParameters,
                                                        &this->_physicalElements);
  }

  void
  assembleOrderParameterForceVector(const vector<ElementVector> & displacements,
                                    const double time,
                                    Eigen::VectorXd * forceVector) const {
    forceVector->fill(0);
    try {
      Assemblers::Utilities::assembleOrderParameterForceVectorUtility(displacements,
                                                                      this->_physicalElements, time,
                                                                      forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the physical elements\n");
      throw;
    }
    try {
      this->_externalForceAssembler.assembleForceVector(displacements, time,
                                                        forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the external force assembler\n");
      throw;
    }
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMobilityMatrix() const {
    const size_t numberOfDofs = this->_numberOfNodes * OrderParameterDofs;
    Eigen::SparseMatrix<double> mobilityMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleLumpedMobilityMatrixUtility(this->_numberOfNodes,
                                                               this->_physicalElements,
                                                               &mobilityMatrix);
    return mobilityMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMobilityMatrix() const {
    const size_t numberOfDofs = this->_numberOfNodes * OrderParameterDofs;
    Eigen::SparseMatrix<double> mobilityMatrix(numberOfDofs, numberOfDofs);
    mobilityMatrix.setZero();
    Assemblers::Utilities::assembleConsistentMobilityMatrixUtility(this->_numberOfNodes,
                                                                   this->_physicalElements,
                                                                   &mobilityMatrix);
    return mobilityMatrix;
  }

};

template <class PhysicalElementType,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType> >
class ElectroMechanicalAssembler : public Assembler<PhysicalElementType, ExternalForceAssembler> {

  typedef Assembler<PhysicalElementType, ExternalForceAssembler> Base;
public:
  typedef typename PhysicalElementType::Vector          ElementVector;
  static const unsigned int MechanicalDegreesOfFreedom = PhysicalElementType::MechanicalDegreesOfFreedom;

  ElectroMechanicalAssembler(const size_t numberOfNodes) : Base(numberOfNodes) {
  }

  template <class QuadratureRule, class MaterialModel>
  ElectroMechanicalAssembler(const SingleElementMesh<PhysicalElementType> & mesh,
                             const typename PhysicalElementType::Properties & elementProperties,
                             const QuadratureRule & quadratureRule,
                             const MaterialModel & materialModel) :
    Base(mesh, elementProperties, quadratureRule, materialModel) {
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMassMatrix() const {
    size_t numberOfDofs = this->_numberOfNodes * MechanicalDegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleMechanicalLumpedMassMatrixUtility(this->_numberOfNodes,
                                                                     this->_physicalElements,
                                                                     &massMatrix);
    return massMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMassMatrix() const {
    size_t numberOfDofs = this->_numberOfNodes * MechanicalDegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs, numberOfDofs);
    massMatrix.setZero();
    Assemblers::Utilities::assembleMechanicalConsistentMassMatrixUtility(this->_numberOfNodes,
                                                                         this->_physicalElements,
                                                                         &massMatrix);
    return massMatrix;
  }

};
*/
	
	
template <class VolumeElement,
class SurfaceElementType1,
class SurfaceElementType2>
class AssemblerChemoMechanicalProblem {
		
	public:
  //typedef typename VolumeElement::MechanicalVector					MechanicalVector;
  //typedef typename VolumeElement::MechanicalStiffnessMatrix MechanicalStiffnessMatrix;
  //typedef typename VolumeElement::ChemicalVector						ChemicalVector;
	typedef typename VolumeElement::PotentialVector						PotentialVector;
	typedef typename VolumeElement::TotalVariableVector				TotalVariableVector;
	//typedef typename VolumeElement::MechanicalForces					MechanicalForces;
	//typedef typename VolumeElement::ChemicalForces						ChemicalForces;
	//typedef typename VolumeElement::PotentialForces						PotentialForces;
	typedef typename VolumeElement::Forces										Forces;
	typedef typename VolumeElement::PrimitiveVariables				PrimitiveVariables;
  typedef typename VolumeElement::Point											ElementPoint;
  typedef typename VolumeElement::Node											ElementNode;
  static const unsigned int VolumeElementNumberOfNodes = VolumeElement::NumberOfNodes;
  static const unsigned int SurfaceElement1NumberOfNodes = SurfaceElementType1::NumberOfNodes;
  static const unsigned int SurfaceElement2NumberOfNodes = SurfaceElementType2::NumberOfNodes;
  static const unsigned int SpatialDimension     = VolumeElement::SpatialDimension;
  static const unsigned int TotalDegreesOfFreedom     = VolumeElement::TotalDegreesOfFreedom;
	//static const unsigned int MechanicalDegreesOfFreedom = VolumeElement::MechanicalDegreesOfFreedom;
		
  AssemblerChemoMechanicalProblem(const size_t numberOfNodes) : _numberOfNodes(numberOfNodes) {
	}
		
  // add all elements
  void
  addElement(const VolumeElement & e) {
		_allVolumeElements.push_back(e);
	}
		
  void
  addElement(const SurfaceElementType1 & e) {
		_allSurfaceType1Elements.push_back(e);
	}
		
  void
  addElement(const SurfaceElementType2 & e) {
		_allSurfaceType2Elements.push_back(e);
	}
	
  void
  assembleForceVector(const vector<TotalVariableVector> & allPrimitives,
											const double phi,
											const double time,
											Eigen::VectorXd * forceVector) const {
		
		Assemblers::Utilities::assembleForceVectorUtility(allPrimitives,phi,
																											_allSurfaceType1Elements,time,
																											forceVector);
		
		Assemblers::Utilities::assembleForceVectorUtility(allPrimitives,phi,
																											_allSurfaceType2Elements,time,
																											forceVector);
		
		Assemblers::Utilities::assembleForceVectorUtility(allPrimitives,phi,
																											_allVolumeElements, time,
																											forceVector);
	}
	
  void
  assembleStiffnessMatrix(const vector<TotalVariableVector> & allPrimitives,
													const double phi,
													const double time,
													Eigen::SparseMatrix<double> * stiffnessMatrix) const {
		Assemblers::Utilities::assembleStiffnessMatrixUtility(allPrimitives,phi,
																													_allVolumeElements,time,
																													stiffnessMatrix);
		Assemblers::Utilities::assembleStiffnessMatrixUtility(allPrimitives,phi,
																													_allSurfaceType1Elements,time,
																													stiffnessMatrix);
		Assemblers::Utilities::assembleStiffnessMatrixUtility(allPrimitives,phi,
																													_allSurfaceType2Elements,time,
																													stiffnessMatrix);
	}
	
	void
	assembleForceVectorOfPreviousConcentrations(const vector<TotalVariableVector> & allPrimitives,
																							const double phi,
																							const double timeStep,
																							Eigen::VectorXd * forceVector) const {
		Assemblers::Utilities::assemblePreviousConcentrationForceVectorUtility(allPrimitives,phi,
																											 _allVolumeElements, timeStep,
																											 forceVector);
	}
	
		
  void
  updateInternalVariables(const vector<TotalVariableVector> & allPrimitives,
													const double time) {
		Assemblers::Utilities::updateInternalVariablesUtility(allPrimitives,time,
																													& _allVolumeElements);
	}

	void
	getAverageConcentration(const vector<TotalVariableVector> & allPrimitives,
													const double timeStep,
													double * averageConcentration) const {
		Assemblers::Utilities::assembleAverageConcentrationUtility(allPrimitives,
																															 _allVolumeElements,
																															 timeStep,
																															 averageConcentration);
	}
	
  const vector<VolumeElement> &
  getPrimaryPhysicalElements() const {
		return _allVolumeElements;
	}
		
  const vector<SurfaceElementType1> &
  getSecondaryPhysicalElements() const {
		return _allSurfaceType1Elements;
	}
		
  const vector<SurfaceElementType2> &
  getTertiaryPhysicalElements() const {
		return _allSurfaceType2Elements;
	}
		
  size_t
  getNumberOfNodes() const {
		return _numberOfNodes;
	}
		
  Eigen::VectorXd
  allocateZeroedGlobalForceVector() const {
		const size_t numberOfDOFs = _numberOfNodes * TotalDegreesOfFreedom;
		Eigen::VectorXd temp(numberOfDOFs);
		temp.fill(0);
		return temp;
	}
		
  Eigen::SparseMatrix<double>
  allocateZeroedStiffnessMatrix() const {
		const size_t numberOfDOFs = _numberOfNodes * TotalDegreesOfFreedom;
		Eigen::SparseMatrix<double> temp(numberOfDOFs, numberOfDOFs);
		return temp;
	}
		
	protected:
  vector<VolumeElement> _allVolumeElements;
  vector<SurfaceElementType1> _allSurfaceType1Elements;
  vector<SurfaceElementType2> _allSurfaceType2Elements;
  size_t _numberOfNodes;
};
	
	
}

#endif  // ASSEMBLER_H
