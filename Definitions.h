// -*- C++ -*-
#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <iostream>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#pragma GCC diagnostic pop
#include <cstdarg>

#define EIGEN_SUPERLU_SUPPORT

// complete list of cell types is here:
//  http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
#include <vtkCellType.h>

#include <ads/utility/ParseOptionsArguments.h>


enum VerboseFlag {Quiet, Verbose};

using std::array;
using std::vector;
using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define errorStatement(s, ...)                                  \
  do {                                                          \
    fprintf (stderr, "(%30s:%40s:%4d) -- " s,                   \
             __FILE__, __func__, __LINE__, ##__VA_ARGS__);      \
    fflush (stderr);                                            \
  } while (0)

#define FILENAMEWITHEXTENSION (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

// from here: http://stackoverflow.com/questions/6420194/how-to-print-both-to-stdout-and-file-in-c
void
tee(FILE *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  va_start(ap, fmt);
  vfprintf(f, fmt, ap);
  va_end(ap);
}

typedef Matrix<double, 1, 1> Vector1d;
typedef Matrix<double, 2, 1> Vector2d;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, 4, 1> Vector4d;
typedef array<array<array<array<double,3>,3>,3>,3> Standard3DTangentMatrix;
typedef array<array<array<array<double,2>,2>,2>,2> Standard2DTangentMatrix;

typedef Matrix<double,2,2> StandardDeformationGradient2D;
typedef Matrix<double,3,3> StandardDeformationGradient3D;


template <class P>
struct NodeWithId {
  NodeWithId(const size_t id,
             const P & position) :
    _id(id), _position(position){
  }
  typedef P Point;
  size_t _id;
  P _position;
  NodeWithId() :
    _id(std::numeric_limits<size_t>::max()), _position() {
  }
};

template <class P>
struct NodeWithIdAndPeriodicGroupNumber {
  NodeWithIdAndPeriodicGroupNumber(const size_t id,
                                   const P & position) :
    _id(id), _position(position),
    _periodicGroupNumber(std::numeric_limits<size_t>::max()){
  }
  typedef P Point;
  size_t _id;
  P _position;
  size_t _periodicGroupNumber;
  NodeWithIdAndPeriodicGroupNumber() :
    _id(std::numeric_limits<size_t>::max()), _position(),
    _periodicGroupNumber(std::numeric_limits<size_t>::max()) {
  }
};

enum NodeLocation {Unknown, Volume, Boundary};

template <class P>
struct NodeWithIdPeriodicGroupNumberAndLocation {
  NodeWithIdPeriodicGroupNumberAndLocation(const size_t id,
                                           const P & position) :
    _id(id), _position(position),
    _periodicGroupNumber(std::numeric_limits<size_t>::max()),
    _nodeLocation(Unknown){
  }
  typedef P Point;
  size_t _id;
  P _position;
  size_t _periodicGroupNumber;
  NodeLocation _nodeLocation;
  NodeWithIdPeriodicGroupNumberAndLocation() :
    _id(std::numeric_limits<size_t>::max()), _position(),
    _periodicGroupNumber(std::numeric_limits<size_t>::max()),
    _nodeLocation(Unknown){
  }
};

template <unsigned D, unsigned QP>
struct QuadratureRule {
  typedef array<Matrix<double, D, 1>, QP> QuadPoints;
  typedef array<double, QP> QuadWeights;
  QuadratureRule(const array<Matrix<double, D, 1>, QP> & points,
                 const array<double, QP> & weights) :
    _points(points), _weights(weights) {
  }
  const array<Matrix<double, D, 1>, QP> _points;
  const array<double, QP> _weights;
};

template <class Element>
struct SingleElementMesh {

  typedef typename Element::Node Node;
  static const unsigned int SpatialDimension = Element::SpatialDimension;

  vector<typename Element::Node> _nodes;
  vector<array<size_t, Element::NumberOfNodes> > _connectivity;
};

template <class Element1, class Element2>
struct TwoElementMesh {

  typedef typename Element1::Node Node;
  static const unsigned int SpatialDimension = Element1::SpatialDimension;

  vector<typename Element1::Node> _nodes;
  vector<array<size_t, Element1::NumberOfNodes> > _connectivity1;
  vector<array<size_t, Element2::NumberOfNodes> > _connectivity2;
};

class EssentialBoundaryCondition {
public:
  size_t _nodeId;
  unsigned int _coordinate;
  double _constraint;
  EssentialBoundaryCondition(const size_t nodeId,
                             const unsigned int coordinate,
                             const double constraint) :
    _nodeId(nodeId), _coordinate(coordinate), _constraint(constraint) {
  }
};

namespace MaterialModels {
class EmptyInternalVariables {
public:
  EmptyInternalVariables(){
  }

  static
  EmptyInternalVariables
  generateRandomAdmissibleTestVariables() {
    return EmptyInternalVariables();
  }
};
}

// this weird looking function is to avoid unused variable warnings from the
//  compiler with -Wall on
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <typename Variable>
void ignoreUnusedVariable(Variable dummy) {
}
template <typename T>
void ignoreUnusedVariables(const T & t) {
}
template <typename T, typename U>
void ignoreUnusedVariables(const T & t, const U & u) {
}
template <typename T, typename U, typename V>
void ignoreUnusedVariables(const T & t, const U & u, const V & v) {
}
template <typename T, typename U, typename V, typename W>
void ignoreUnusedVariables(const T & t, const U & u, const V & v, const W & w) {
}
#pragma GCC diagnostic pop

struct SparseEigenMatrixRowZeroer {
  // constructor
  SparseEigenMatrixRowZeroer(const vector<size_t> & rowIndices) {
    std::copy(rowIndices.begin(), rowIndices.end(),
              std::inserter(_rowIndices, _rowIndices.begin()));
  }
  bool operator() (const size_t & row, const size_t & col,
                   const double & value) const{
    ignoreUnusedVariables(col, value);
    return _rowIndices.find(row) == _rowIndices.end();
  };
private:
  std::set<size_t> _rowIndices;
};

struct SparseEigenMatrixColZeroer {
  // constructor
  SparseEigenMatrixColZeroer(const vector<size_t> & colIndices) {
    std::copy(colIndices.begin(), colIndices.end(),
              std::inserter(_colIndices, _colIndices.begin()));
  }
  bool operator() (const size_t & row, const size_t & col,
                   const double & value) const{
    ignoreUnusedVariables(row, value);
    return _colIndices.find(col) == _colIndices.end();
  };
private:
  std::set<size_t> _colIndices;
};

struct SparseEigenMatrixRemoveNumbersUnderAbsThreashold{
  SparseEigenMatrixRemoveNumbersUnderAbsThreashold(const double lowerAbsThreashold) :
    _lowerAbsThreashold(lowerAbsThreashold){
  }

  bool operator() (const size_t & row, const size_t & col,
                   const std::complex<double> & value) const{
    ignoreUnusedVariables(row,col);
    const bool keepValue = (std::abs(value) < _lowerAbsThreashold ? false: true);
    return keepValue;
  }
private:
  double _lowerAbsThreashold;
};

char exceptionBuffer[10000];
#define throwException(s, ...)                  \
  sprintf(exceptionBuffer, s, ##__VA_ARGS__);   \
  throw std::runtime_error(exceptionBuffer);

double
smoothHeaviside(const double x,
                const double heavisideTransitionPoint,
                const double sharpnessOfTransition = 1){
  return 1./double(1+exp(-2*sharpnessOfTransition * (x - heavisideTransitionPoint)));
}

template<class Node, size_t NumberOfNodesPerBoundary>
struct SingleElementBoundary{
  SingleElementBoundary(array<Node,NumberOfNodesPerBoundary> nodes):
    _nodes(nodes),
    _centerOfBoundary(){
    _centerOfBoundary.setZero();
    for (size_t nodeIndex = 0; nodeIndex < NumberOfNodesPerBoundary; nodeIndex++){
      _centerOfBoundary += _nodes[nodeIndex]._position;
    }
    _centerOfBoundary /= double(NumberOfNodesPerBoundary);
  }
  array<Node,NumberOfNodesPerBoundary> _nodes;
  typename Node::Point _centerOfBoundary;

  SingleElementBoundary() :
    _nodes(),
    _centerOfBoundary(){
    _centerOfBoundary.setZero();
  }
};

template<class SingleElementBoundary,size_t NumberOfBoundaries>
struct AllElementBoundaries{
  array<SingleElementBoundary,NumberOfBoundaries> _boundaries;
};

#endif // DEFINITIONS_H
