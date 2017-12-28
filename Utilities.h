// -*- C++ -*-
#ifndef UTILITIES_H
#define UTILITIES_H

#include "Definitions.h"
#include <Eigen/Eigenvalues>
#include <sys/stat.h>
#include <mpi.h>


namespace Utilities {

template <class Element>
array<typename Element::Vector, Element::NumberOfNodes>
getElementDisplacementsFromGlobalList(const array<size_t, Element::NumberOfNodes> & elementNodeIds,
                                      const vector<typename Element::Vector> & displacements) {
  array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements;
  for (size_t nodeIndex = 0;
       nodeIndex < Element::NumberOfNodes; ++nodeIndex) {
    size_t nodeId = elementNodeIds[nodeIndex];
    // check that the nodeId is valid
    if (nodeId >= displacements.size() || nodeId < 0) {
      errorStatement("Utilities::getElementDisplacementsFromGlobalList found an error "
                     "found a node Id of %zu but there are only %zu "
                     "provided displacements\n", nodeId, displacements.size());
      // TODO: throw an exception
      assert(nodeId < displacements.size() && nodeId > 0);
    }
    elementDisplacements[nodeIndex] = displacements[nodeId];
  }
  return elementDisplacements;
}
	
template <class Element>
array<typename Element::TotalVariableVector, Element::NumberOfNodes>
getElementPrimitivesFromGlobalList(const array<size_t, Element::NumberOfNodes> & elementNodeIds,
																		 const vector<typename Element::TotalVariableVector> & primitives) {
  array<typename Element::TotalVariableVector, Element::NumberOfNodes> elementPrimitives;
  for (size_t nodeIndex = 0;
			 nodeIndex < Element::NumberOfNodes; ++nodeIndex) {
		size_t nodeId = elementNodeIds[nodeIndex];
		// check that the nodeId is valid
		if (nodeId >= primitives.size() || nodeId < 0) {
			errorStatement("Utilities::getElementDisplacementsFromGlobalList found an error "
										 "found a node Id of %zu but there are only %zu "
										 "provided displacements\n", nodeId, primitives.size());
			// TODO: throw an exception
			assert(nodeId < primitives.size() && nodeId > 0);
		}
		elementPrimitives[nodeIndex] = primitives[nodeId];
	}
  return elementPrimitives;
}
	

template <class Element>
array<typename Element::OrderParameter, Element::NumberOfNodes>
getElementOrderParametersFromGlobalList(const array<size_t, Element::NumberOfNodes> & elementNodeIds,
                                        const vector<typename Element::OrderParameter> & orderParameters) {
  array<typename Element::OrderParameter, Element::NumberOfNodes> elementOrderParameters;
  for (size_t nodeIndex = 0;
       nodeIndex < Element::NumberOfNodes; ++nodeIndex) {
    size_t nodeId = elementNodeIds[nodeIndex];
    // check that the nodeId is valid
    if (nodeId >= orderParameters.size() || nodeId < 0) {
      errorStatement("Utilities::getElementDisplacementsFromGlobalList found an error "
                     "found a node Id of %zu but there are only %zu "
                     "provided displacements\n", nodeId, orderParameters.size());
      // TODO: throw an exception
      assert(nodeId < orderParameters.size() && nodeId > 0);
    }
    elementOrderParameters[nodeIndex] = orderParameters[nodeId];
  }
  return elementOrderParameters;
}

template <class Assembler>
vector<typename Assembler::ElementVector>
distributeGlobalVectorToLocalVectors(const Eigen::VectorXd & globalVector) {
  const unsigned int DegreesOfFreedom = Assembler::DegreesOfFreedom;

  unsigned int numberOfNodes = globalVector.size() / DegreesOfFreedom;
  if (numberOfNodes*DegreesOfFreedom != globalVector.size()) {
    printf("Cannot split a global vector of size %zu into local vectors "
           "for elements of spatial dimension %u\n", globalVector.size(),
           DegreesOfFreedom);
  }

  vector<typename Assembler::ElementVector> localVectors(numberOfNodes);
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
    for (size_t coordinate = 0; coordinate < DegreesOfFreedom; coordinate++) {
      localVectors[nodeIndex](coordinate) =
        globalVector(nodeIndex * DegreesOfFreedom + coordinate);
    }
  }
  return localVectors;
}
	
template <class Assembler>
vector<typename Assembler::TotalVariableVector>
distributeGlobalVectorToLocalVectors(const Eigen::VectorXd & globalVector) {
  const unsigned int DegreesOfFreedom = Assembler::TotalDegreesOfFreedom;
		
  unsigned int numberOfNodes = globalVector.size() / DegreesOfFreedom;
  if (numberOfNodes*DegreesOfFreedom != globalVector.size()) {
		printf("Cannot split a global vector of size %zu into local vectors "
					 "for elements of spatial dimension %u\n", globalVector.size(),
					 DegreesOfFreedom);
	}
		
  vector<typename Assembler::TotalVariableVector> localVectors(numberOfNodes);
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
		for (size_t coordinate = 0; coordinate < DegreesOfFreedom; coordinate++) {
			localVectors[nodeIndex](coordinate) =
			globalVector(nodeIndex * DegreesOfFreedom + coordinate);
		}
	}
  return localVectors;
}

template <class Assembler>
vector<typename Assembler::OrderParameter>
distributeGlobalOrderParameterVectorToLocalVectors(const Eigen::VectorXd & globalVector) {
  const unsigned int DegreesOfFreedom = Assembler::OrderParameterDofs;

  unsigned int numberOfNodes = globalVector.size() / DegreesOfFreedom;
  if (numberOfNodes*DegreesOfFreedom != globalVector.size()) {
    printf("Cannot split a global vector of size %zu into local vectors "
           "for elements of spatial dimension %u\n", globalVector.size(),
           DegreesOfFreedom);
  }

  vector<typename Assembler::OrderParameter> localVectors(numberOfNodes);
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
    for (size_t coordinate = 0; coordinate < DegreesOfFreedom; coordinate++) {
      localVectors[nodeIndex](coordinate) =
        globalVector(nodeIndex * DegreesOfFreedom + coordinate);
    }
  }
  return localVectors;
}

void
computeEigenSolution(const MatrixXd & A, const MatrixXd & B,
                     vector<double> * eigenvalues, vector<VectorXd> * eigenvectors) {

  if (A.rows() != B.rows() || A.cols() != B.cols() || A.rows() != A.cols()) {
    printf("Square matrix arguments A and B to computeEigenSolution must have "
           "the same size, but A is %zdx%zd and B is %zdx%zd\n",
           A.rows(), A.cols(), B.rows(), B.cols());
    exit(1);
  }

  Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> es(A, B, Eigen::ComputeEigenvectors | Eigen::Ax_lBx);

  VectorXd eigenvaluesVector = es.eigenvalues();
  MatrixXd eigenvectorsMatrix = es.eigenvectors();
  eigenvalues->resize(A.cols());
  eigenvectors->resize(A.cols());
  for (size_t colIndex = 0; colIndex < size_t( A.cols()); ++colIndex) {
    (*eigenvalues)[colIndex] = eigenvaluesVector(colIndex);
    (*eigenvectors)[colIndex] = eigenvectorsMatrix.col(colIndex);
  }
}

template <class Element>
array<typename Element::Vector,Element::NumberOfNodes>
packArrays(const Matrix<double, (Element::NumberOfNodes) * (Element::DegreesOfFreedom), 1> & unpacked) {
  array<typename Element::Vector,Element::NumberOfNodes> packed;
  for (unsigned int nodeIndex = 0; nodeIndex < Element::NumberOfNodes; nodeIndex++) {
    for (unsigned int dofIndex = 0; dofIndex < Element::DegreesOfFreedom; dofIndex++) {
      packed[nodeIndex](dofIndex) = unpacked(Element::DegreesOfFreedom*nodeIndex+dofIndex);
    }
  }
  return packed;
}

template <class Element>
array<typename Element::OrderParameter,Element::NumberOfNodes>
packArraysOfOrderParameters(const Matrix<double, (Element::NumberOfNodes) * (Element::OrderParameterDofs), 1> & unpacked) {
  array<typename Element::OrderParameter,Element::NumberOfNodes> packed;
  for (unsigned int nodeIndex = 0; nodeIndex < Element::NumberOfNodes; nodeIndex++) {
    for (unsigned int dofIndex = 0; dofIndex < Element::OrderParameterDofs; dofIndex++) {
      packed[nodeIndex](dofIndex) = unpacked(Element::OrderParameterDofs*nodeIndex+dofIndex);
    }
  }
  return packed;
}

// Function for extracting out particular degrees of freedom in an element from a global list.
// The ordering of degrees of freedom is the same as the order they appear in degreeOfFreedomIndices.
template <class Element, unsigned NumberOfDegreesOfFreedomToExtract>
void extractDegreesOfFreedomFromGlobalList(const vector<typename Element::Vector> & displacements,
                                           const array<unsigned int, NumberOfDegreesOfFreedomToExtract> & degreeOfFreedomIndices,
                                           VectorXd * extractedDegreesOfFreedom) {

  const size_t numberOfElements = displacements.size();

  if (NumberOfDegreesOfFreedomToExtract == 0) {
    cout << "Error: no indices provided." << endl;
    exit(-1);
  }
  if (extractedDegreesOfFreedom->size() != numberOfElements*NumberOfDegreesOfFreedomToExtract) {
    cout << "Error: output vector wrong size." << endl;
    exit(-1);
  }

  for (unsigned int elementIndex=0; elementIndex<numberOfElements; elementIndex++) {
    for (unsigned int dofIndex=0; dofIndex<NumberOfDegreesOfFreedomToExtract; dofIndex++) {
      const unsigned int degreeOfFreedom = degreeOfFreedomIndices[dofIndex];
      if (degreeOfFreedom > Element::DegreesOfFreedom) {
        cout << "Error: invalid index." << endl;
        exit(-1);
      }
      (*extractedDegreesOfFreedom)(NumberOfDegreesOfFreedomToExtract*elementIndex + dofIndex)
        = displacements[elementIndex](degreeOfFreedom);
    }
  }

}

//**********Code below is from before the great code shift of Spring 2013*********//


template <unsigned Nodes, unsigned Dofs>
inline
Matrix<double, Nodes*Dofs, 1>
unpackArrays(const array<Matrix<double, Dofs, 1>, Nodes> & arrays) {
  Matrix<double, Nodes*Dofs, 1> unpacked;
  for (size_t i = 0; i < Nodes; ++i) {
    for (size_t j = 0; j < Dofs; ++j) {
      unpacked(i*Dofs+j) = arrays[i](j);
    }
  }
  return unpacked;
}

template <class Element>
inline
Matrix<double, Element::NumberOfNodes*Element::DegreesOfFreedom, 1>
unpackArrays(const array<Matrix<double, Element::DegreesOfFreedom, 1>, Element::NumberOfNodes> & arrays) {
  return unpackArrays<Element::NumberOfNodes, Element::DegreesOfFreedom>(arrays);
}

template <unsigned Nodes, unsigned Dofs>
inline
array<Matrix<double, Dofs, 1>, Nodes>
packArrays(const Matrix<double, Nodes*Dofs, 1> & unpacked) {
  array<Matrix<double, Dofs, 1>, Nodes> arrays;
  for (size_t i = 0; i < Nodes; ++i) {
    for (size_t j = 0; j < Dofs; ++j) {
      arrays[i](j) = unpacked(i*Dofs+j);
    }
  }
  return arrays;
}
/*
  template <class Element>
  inline
  array<Matrix<double, Element::DegreesOfFreedom, 1>, Element::NumberOfNodes>
  packArrays(const Matrix<double, Element::NumberOfNodes*Element::DegreesOfFreedom, 1> & unpacked) {
  return packArrays<Element::NumberOfNodes, Element::DegreesOfFreedom>(unpacked);
  }
*/
template <class Element>
void
printStresses(const vector<Element> & elements,
              const vector<typename Element::Vector> & displacements) {
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements =
      Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                displacements);
    const array<typename Element::Stress, Element::QuadPoints> qpStresses =
      element.computeStressesAtGaussPoints(elementDisplacements);
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      printf("%03zu:%02zu - (", elementIndex, qpIndex);
      for (int i = 0; i < qpStresses[qpIndex].size(); ++i) {
        printf("%10.6f", qpStresses[qpIndex](i));
        if (i != qpStresses[qpIndex].size() - 1) {
          printf(", ");
        } else {
          printf(")\n");
        }
      }
    }
  }
}

template <class Element>
typename Element::Stress
computeMaxStressDifferenceFromFirstElement(const vector<Element> & elements,
                                           const vector<typename Element::Vector> & displacements,
                                           const double time) {
  typedef typename Element::Stress Stress;
  Stress correctStress;
  Stress maxDifferences = Stress::Zero();
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements =
      Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                displacements);
    const array<Stress, Element::QuadPoints> qpStresses =
      elements[elementIndex].computeStressesAtGaussPoints(elementDisplacements,
                                                          time);
    if (elementIndex == 0) {
      correctStress = qpStresses[0];
    }
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      for (unsigned int component = 0; component < qpStresses[qpIndex].size();
           ++component) {
        const double thisComponentsDifference =
          std::abs(qpStresses[qpIndex](component) - correctStress(component));
        maxDifferences(component) =
          std::max(thisComponentsDifference,
                   maxDifferences(component));
      }
    }
  }
  return maxDifferences;
}

Eigen::Matrix3d
createRotationMatrix(const Eigen::Vector3d & e, const double theta) {
  Eigen::Vector3d eNormalized;
  eNormalized << e[0], e[1], e[2];
  eNormalized.normalize();
  Eigen::Matrix3d temp;
  temp << 0 , -eNormalized[2] , eNormalized[1], eNormalized[2] , 0 , -eNormalized[0], -eNormalized[1] , eNormalized[0] , 0;
  return (eNormalized*eNormalized.transpose())*(1-cos(theta)) +
    Eigen::Matrix3d::Identity()*cos(theta) + temp*sin(theta);
}

enum TrimStyle {DontTrim, Trim};

// retrieved from http://oxaric.wordpress.com/2008/11/23/3-simple-c-functions/
std::string trim( std::string line ) {
  if ( line.empty() ) {
    return "";
  }
  int string_size = (int)(line.length());
  int beginning_of_string = 0;
  // the minus 1 is needed to start at the first character
  // and skip the string delimiter
  int end_of_string = string_size - 1;
  bool encountered_characters = false;
  // find the start of chracters in the string
  while ( (beginning_of_string < string_size) && (!encountered_characters) ) {
    // if a space or tab was found then ignore it
    if ( (line[ beginning_of_string ] != ' ') &&
         (line[ beginning_of_string ] != '\t') ) {
      encountered_characters = true;
    } else {
      beginning_of_string++;
    }
  }
  // test if no characters were found in the string
  if ( beginning_of_string == string_size ) {
    return "";
  }
  encountered_characters = false;
  // find the character in the string
  while ( (end_of_string > beginning_of_string) &&
          (!encountered_characters) ) {
    // if a space or tab was found then ignore it
    if ( (line[ end_of_string ] != ' ') && (line[ end_of_string ] != '\t') ) {
      encountered_characters = true;
    } else {
      end_of_string--;
    }
  }
  // return the original string with all whitespace removed from
  // its beginning and end
  // + 1 at the end to add the space for the string delimiter
  return line.substr( beginning_of_string,
                      end_of_string - beginning_of_string + 1 );
}

void tokenize(const std::string& str, const std::string& delimiters,
              const TrimStyle doTrim, std::vector<std::string>* tokens) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    std::string candidate = str.substr(lastPos, pos - lastPos);
    if (doTrim == Trim) {
      candidate = trim(candidate);
    }
    tokens->push_back(candidate);
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

template <class T>
T
convertString(const string & str) {
  std::stringstream ss;
  ss << str;
  T returnVal;
  ss >> returnVal;
  return returnVal;
}

template <>
float
convertString<float>(const string & str) {
  if (strcmp(str.c_str(), "nan") == 0) {
    //printf("Parsed a nan\n");
    return std::numeric_limits<float>::quiet_NaN();
  } else {
    std::stringstream ss;
    ss << str;
    float returnVal;
    ss >> returnVal;
    return returnVal;
  }
}

template <>
double
convertString<double>(const string & str) {
  if (strcmp(str.c_str(), "nan") == 0) {
    //printf("Parsed a nan\n");
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    std::stringstream ss;
    ss << str;
    double returnVal;
    ss >> returnVal;
    return returnVal;
  }
}

vector<std::array<unsigned int, 3> >
triangulate2dPoints(const vector<Matrix<double, 2, 1> > & points,
                    const string & triangleExecutablePath) {
  FILE* nodesFile = fopen("nodes.node", "w");
  fprintf(nodesFile, "%zu 2 0 0\n", points.size());
  for (unsigned int index = 0; index < points.size(); ++index) {
    fprintf(nodesFile, "%u %lf %lf\n", index,
            points[index](0),
            points[index](1));
  }
  fclose(nodesFile);

  // run the triangulation
  char buffer[1000];
  sprintf(buffer, "%s nodes.node", triangleExecutablePath.c_str());
  printf("attempting to run command %s\n", buffer);
  system(buffer);

  // parse output
  ifstream elementsFile("nodes.1.ele");
  string line;
  getline(elementsFile, line);
  vector<string> tokens;
  tokenize(line, " ", Trim, &tokens);
  if (tokens.size() != 3) {
    fprintf(stderr, "could not parse triangle ele header, found %zu "
            "tokens in line \"%s\"\n", tokens.size(), line.c_str());
    exit(1);
  }
  const unsigned int numberOfTriangles =
    convertString<unsigned int>(tokens[0]);
  vector<std::array<unsigned int, 3> > triangles;
  for (unsigned int triangleIndex = 0;
       triangleIndex < numberOfTriangles; ++triangleIndex) {
    getline(elementsFile, line);
    tokens.resize(0);
    tokenize(line, " ", Trim, &tokens);
    if (tokens.size() != 4) {
      fprintf(stderr, "could not parse triangle element number %u, found %zu "
              "tokens in line \"%s\"\n", triangleIndex, tokens.size(),
              line.c_str());
      exit(1);
    }
    std::array<unsigned int, 3> triangle;
    triangle[0] = convertString<unsigned int>(tokens[1]);
    triangle[1] = convertString<unsigned int>(tokens[2]);
    triangle[2] = convertString<unsigned int>(tokens[3]);
    triangles.push_back(triangle);
  }
  elementsFile.close();

  return triangles;
}

void
directoryCreator(const string & fullPath,
                 const bool createNewDirectory = true,
                 const VerboseFlag verboseFlag = Verbose){

  // split the fullPath by slashes
  vector<string> tokens;
  Utilities::tokenize(fullPath, "/", Utilities::Trim, &tokens);

  string incrementalPath;
  for (unsigned int tokenIndex = 0; tokenIndex < tokens.size(); ++tokenIndex) {

    // keep building the incrementalPath
    if (tokenIndex > 0) {
      incrementalPath += "/";
    }
    incrementalPath += tokens[tokenIndex];

    std::ifstream ifile(incrementalPath.c_str());
    if ((bool)ifile == false && createNewDirectory == true) {
      if (verboseFlag == Verbose){
        printf("Did not find a directory at %s, creating one.\n",
               incrementalPath.c_str());
      }
      int status =
        mkdir(incrementalPath.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (status != 0 && verboseFlag == Verbose) {
        fprintf(stderr, "Problem creating directory %s\n",
                incrementalPath.c_str());
      }
    } else if ((bool)ifile == false && createNewDirectory == false
               && verboseFlag == Verbose){
      printf("Couldn't find the directory at %s that you specified,"
             "and you said not to make it\n",incrementalPath.c_str());
    } else if (verboseFlag == Verbose){
      printf("Found a directory already existing at %s, not creating one.\n",
             incrementalPath.c_str());
    }
  }
}

template<class MatrixType>
void
checkForZeroRowsOfSparseMatrix(const MatrixType & sparseMatrix,
                               const string matrixName){
  vector<double> collapseColumnsToOne;
  collapseColumnsToOne.reserve(sparseMatrix.rows());
  for (size_t rowIndex = 0; rowIndex < sparseMatrix.rows(); rowIndex++){
    collapseColumnsToOne.push_back(0);
  }
  for (size_t colIndex = 0; colIndex < sparseMatrix.outerSize(); ++colIndex){
    for (typename MatrixType::InnerIterator
           iterator(sparseMatrix,colIndex);
         iterator; ++iterator){
      collapseColumnsToOne[iterator.row()] += std::abs(iterator.value());
    }
  }
  vector<size_t> zeroRows;
  for (size_t rowIndex = 0; rowIndex < sparseMatrix.rows(); rowIndex++){
    if (collapseColumnsToOne[rowIndex] == 0){
      zeroRows.push_back(rowIndex);
    }
  }

  if (zeroRows.size() > 0) {
    printf("Error: %s has (%lu) zero rows:\n", matrixName.c_str(),zeroRows.size());
    for (size_t i = 0; i < zeroRows.size(); ++i) {
      printf("%zu, ", zeroRows[i]);
    }
    throwException("there were rows of zeros in the checked sparse matrix\n");
  }
}

template<class MatrixType>
void
checkForZeroColsOfSparseMatrix(const MatrixType & sparseMatrix,
                               const string matrixName){
  vector<double> collapseRowsToOne;
  collapseRowsToOne.reserve(sparseMatrix.cols());
  for (size_t colIndex = 0; colIndex < sparseMatrix.cols(); colIndex++){
    collapseRowsToOne.push_back(0);
  }
  for (size_t colIndex = 0; colIndex < sparseMatrix.outerSize(); ++colIndex){
    for (typename MatrixType::InnerIterator
           iterator(sparseMatrix,colIndex);
         iterator; ++iterator){
      collapseRowsToOne[iterator.col()] += std::abs(iterator.value());
    }
  }
  vector<size_t> zeroCols;
  for (size_t colIndex = 0; colIndex < sparseMatrix.cols(); colIndex++){
    if (collapseRowsToOne[colIndex] == 0){
      zeroCols.push_back(colIndex);
    }
  }

  if (zeroCols.size() > 0) {
    printf("Error: %s has (%lu) zero cols:\n", matrixName.c_str(),zeroCols.size());
    for (size_t i = 0; i < zeroCols.size(); ++i) {
      printf("%zu, ", zeroCols[i]);
    }
    throwException("there were cols of zeros in the checked sparse matrix\n");
  }
}

string
getLocalTimeString() {
  char timeBuffer[100];
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  strftime (timeBuffer,80,"%X",timeinfo);
  return string(timeBuffer);
}

void
initializeMpi(int* argc, char*** argv, unsigned int * rank,
              unsigned int * numberOfProcesses) {
  int mpiErrorCode = MPI_Init(argc, argv);
  if (mpiErrorCode != MPI_SUCCESS) {
    char errorString[1000];
    int errorStringLength;

    MPI_Error_string(mpiErrorCode, errorString, &errorStringLength);
    printf("There was an mpi error: %s", errorString);
    exit(1);
  }

  int temp;
  MPI_Comm_rank(MPI_COMM_WORLD, &temp);
  *rank = temp;
  MPI_Comm_size(MPI_COMM_WORLD, &temp);
  *numberOfProcesses = temp;
}

void
initializeMpi(int* argc, char*** argv) {
  unsigned int rank;
  unsigned int numberOfProcesses;
  initializeMpi(argc, argv, &rank, &numberOfProcesses);
}

// A function for computing the apparant stiffness matrix of a mesh for a
// given set of applied nodal forces and displacement boundary conditions.
// Zero body forces are assumed.
template <class Assembler>
Eigen::MatrixXcd
computeEffectiveStiffnessMatrix(vector<typename Assembler::ElementVector> & nodalDisplacements,
    const double time, const vector<EssentialBoundaryCondition> & essentialBCs, 
    const vector<EssentialBoundaryCondition> & forceBCs, Assembler * assembler) {

  const unsigned int numberOfDOFs = assembler->getNumberOfNodes() * Assembler::DegreesOfFreedom;
  const size_t numberOfAppliedForces = forceBCs.size();
  Eigen::SparseMatrix<std::complex<double> > tangentMatrix(numberOfDOFs, numberOfDOFs);

  // have the assembler compute the stiffness matrix
  assembler->assembleStiffnessMatrix(nodalDisplacements, time,  &tangentMatrix);
  Eigen::MatrixXcd fullTangentMatrix = tangentMatrix;

  // apply essential dispalcement boundary conditions to stiffness matrix
  // (zero out row and then put a 1 on the diagonal)
  for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
    const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
    const size_t dofIndex = bc._nodeId * Assembler::DegreesOfFreedom + bc._coordinate;
    for (size_t columnIndex = 0; columnIndex < numberOfDOFs; ++columnIndex) {
      fullTangentMatrix(dofIndex,columnIndex) = 0.;
    }
    fullTangentMatrix(dofIndex,dofIndex) = 1.;
  }


  // Write matrix in block form:
  // [ Kaa Kab ][Ua]   [Fa]
  // [ Kab Kbb ][Ub] = [ 0],
  //
  // where the right hand side has the applied nodal forces on the top
  // and zero external forces are applied on the rest of the nodes.

  // Make Kaa matrix:
  Eigen::MatrixXcd upperLeftBlock(numberOfAppliedForces,numberOfAppliedForces);
  upperLeftBlock.fill(0.);

  size_t blockIndex_i = 0;
  size_t blockIndex_j = 0;
  for (size_t forceBCIndex_i = 0; forceBCIndex_i < numberOfAppliedForces; forceBCIndex_i++) {
    const EssentialBoundaryCondition & bc_i = forceBCs[forceBCIndex_i];
    const size_t dofIndex_i = bc_i._nodeId * Assembler::DegreesOfFreedom + bc_i._coordinate;
    //force_i = bc_i._value;;
    for (size_t forceBCIndex_j = 0; forceBCIndex_j < numberOfAppliedForces; forceBCIndex_j++) {
      const EssentialBoundaryCondition & bc_j = forceBCs[forceBCIndex_j];
      const size_t dofIndex_j = bc_j._nodeId * Assembler::DegreesOfFreedom + bc_j._coordinate;
      //force_j = bc_j._value;;
      upperLeftBlock(blockIndex_i,blockIndex_j) = fullTangentMatrix(dofIndex_i,dofIndex_j);
      blockIndex_j++;
    }
    blockIndex_i++;
    blockIndex_j = 0;
  }

  // Make Kab matrix (Kba is the transpose):
  Eigen::MatrixXcd upperRightBlock(numberOfAppliedForces,numberOfDOFs - numberOfAppliedForces);
  upperRightBlock.fill(0.);

  blockIndex_i = 0;
  blockIndex_j = 0;
  for (size_t forceBCIndex_i = 0; forceBCIndex_i < numberOfAppliedForces; forceBCIndex_i++) {
    const EssentialBoundaryCondition & bc_i = forceBCs[forceBCIndex_i];
    const size_t dofIndex_i = bc_i._nodeId * Assembler::DegreesOfFreedom + bc_i._coordinate;
    for (size_t freeDofIndex = 0; freeDofIndex < numberOfDOFs; freeDofIndex++) {
      bool isFreeIndex = 1;
      for (size_t forceBCIndex_j = 0; forceBCIndex_j < numberOfAppliedForces; forceBCIndex_j++) {
        const EssentialBoundaryCondition & bc_j = forceBCs[forceBCIndex_j];
        const size_t dofIndex_j = bc_j._nodeId * Assembler::DegreesOfFreedom + bc_j._coordinate;
        // check if freeDofIndex is an applied force index
        isFreeIndex = isFreeIndex && (!(freeDofIndex == dofIndex_j));
      }
      if (isFreeIndex) {
        upperRightBlock(blockIndex_i,blockIndex_j) = fullTangentMatrix(dofIndex_i,freeDofIndex);
        blockIndex_j++;
      }
    }
    blockIndex_i++;
    blockIndex_j = 0;
  }

  // Make Kbb matrix:
  Eigen::MatrixXcd bottomRightBlock(numberOfDOFs - numberOfAppliedForces,numberOfDOFs - numberOfAppliedForces);
  bottomRightBlock.fill(0.);

  blockIndex_i = 0;
  blockIndex_j = 0;
  for (size_t freeDofIndex_i = 0; freeDofIndex_i < numberOfDOFs; freeDofIndex_i++) {
    bool isFreeIndex_i = 1;
    for (size_t forceBCIndex_i = 0; forceBCIndex_i < numberOfAppliedForces; forceBCIndex_i++) {
      const EssentialBoundaryCondition & bc_i = forceBCs[forceBCIndex_i];
      const size_t dofIndex_i = bc_i._nodeId * Assembler::DegreesOfFreedom + bc_i._coordinate;
      // check if freeDofIndex is an applied force index
      isFreeIndex_i = isFreeIndex_i && (!(freeDofIndex_i == dofIndex_i));
    }
    if (isFreeIndex_i) {
      for (size_t freeDofIndex_j = 0; freeDofIndex_j < numberOfDOFs; freeDofIndex_j++) {
        bool isFreeIndex_j = 1;
        for (size_t forceBCIndex_j = 0; forceBCIndex_j < numberOfAppliedForces; forceBCIndex_j++) {
          const EssentialBoundaryCondition & bc_j = forceBCs[forceBCIndex_j];
          const size_t dofIndex_j = bc_j._nodeId * Assembler::DegreesOfFreedom + bc_j._coordinate;
          // check if freeDofIndex is an applied force index
          isFreeIndex_j = isFreeIndex_j && (!(freeDofIndex_j == dofIndex_j));
        }
        if (isFreeIndex_j) {
          bottomRightBlock(blockIndex_i,blockIndex_j) = fullTangentMatrix(freeDofIndex_i,freeDofIndex_j);
          blockIndex_j++;
        }
      }
      blockIndex_i++;
      blockIndex_j = 0;
    }
  }

  // We wish to write the expression for:
  // [Keff][Ua] = [Fa],
  // where
  // Keff = Kaa - Kab*(Kbb^-1)*Kba;

  // Compute Kba by transposing Kab
  Eigen::MatrixXcd bottomLeftBlock = upperRightBlock.transpose();

  // Compute (Kbb^-1)*Kba
  Eigen::MatrixXcd sol = bottomRightBlock.fullPivLu().solve(bottomLeftBlock);

  // Compute Keff
  Eigen::MatrixXcd effectiveStiffnessMatrix = upperLeftBlock - upperRightBlock*sol;

  return effectiveStiffnessMatrix;

}

}

#endif  // UTILITIES_H
