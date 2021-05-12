#include <algorithm>
#include <complex>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <Eigen/Dense>
#include <boost/container/static_vector.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/mp11/algorithm.hpp>
using namespace std;
using namespace Eigen;
using namespace boost::mp11;
using boost::container::static_vector;

#ifndef MIN_VERTS
#define MIN_VERTS 2
#endif

#ifndef MAX_VERTS
#define MAX_VERTS 6
#endif

#ifndef VALUE
#define VALUE float
#endif

#ifndef MAX_ITERS
#define MAX_ITERS 500
#endif

typedef VALUE Value;
typedef complex<Value> Complex;

#ifndef EPS1
#define EPS1 (is_same<float, Value>::value ? 1e-6 : 1e-10)
#endif

#ifndef EPS2
#define EPS2 (is_same<float, Value>::value ? 1e-3 : 1e-8)
#endif

static string line;
static ostream* errStream;
static unsigned long numErrs;

static constexpr Value slope(const Complex& v) {
  return v.imag() / v.real();
}

template <typename HullVector, typename PointVector>
void cHullSym(HullVector& modHull, const PointVector& points, Value eps) {
  static_vector<Complex, MAX_VERTS> sPoints;
  static_vector<Value, MAX_VERTS> slopes;
  for (auto it = points.cbegin(); it != points.cend(); it++)
    if (it->imag() > -eps) {
      auto it2 = sPoints.begin();
      while (it2 != sPoints.end()) {
#ifdef SORT_EPS
        if (it2->real() > it->real() - eps) {
          if (it2->real() >= it->real() + eps) { /* braces hack */ <%}<%}
#else
        if (it2->real() >= it->real()) {
          if (it2->real() > it->real()) {
#endif
            sPoints.insert(it2, *it);
          } else if (it->imag() > it2->imag()) {
            it2->imag(it->imag());
          }
          break;
        }
        it2++;
      }
      if (it2 == sPoints.cend())
        sPoints.push_back(*it);
    }
  modHull.clear();
  modHull.push_back(sPoints[0]);
  if (sPoints.size() == 1)
    return;
  modHull.push_back(sPoints[1]);
  slopes.push_back(numeric_limits<Value>::infinity());
  slopes.push_back(slope(sPoints[1] - sPoints[0]));
  int p = 2;
  while (p < sPoints.size()) {
    Value cSlope = slope(sPoints[p] - modHull.back());
    if (cSlope <= slopes.back() + eps) {
      modHull.push_back(sPoints[p]);
      slopes.push_back(cSlope);
      p++;
    } else {
      modHull.pop_back();
      slopes.pop_back();
    }
  }
}

template <int MSize>
bool nrIsPolygon(const Matrix<Value, MSize, MSize>& iMat, Value eps1, Value eps2) {
  typedef Matrix<Complex, MSize, MSize> CXMatrix;
  typedef Matrix<Value, MSize, MSize> VXMatrix;
  EigenSolver<VXMatrix> vEigen;
  static_vector<Complex, MAX_VERTS> hull;
  CXMatrix testMatrix, mulMatrix, cIMat;
  VXMatrix vTestMatrix;
  vEigen.setMaxIterations(MAX_ITERS);
  vEigen.compute(iMat, false);
  if (vEigen.info() != Success) {
    (*errStream) << line << endl;
    numErrs++;
    return false;
  }
  auto eVals_init = vEigen.eigenvalues();
  cHullSym(hull, eVals_init, eps1);
  vTestMatrix = (iMat + iMat.adjoint()) / 2;
  LLT<VXMatrix> loLLT(vTestMatrix - (hull[0].real() - eps2) * VXMatrix::Identity());
  if (loLLT.info() != Success)
    return false;
  LLT<VXMatrix> hiLLT((hull.back().real() + eps2) * VXMatrix::Identity() - vTestMatrix);
  if (hiLLT.info() != Success)
    return false;
  for (int i = 0; i < MSize; i++)
    for (int j = 0; j < MSize; j++)
      cIMat(i,j).real(iMat(i,j));
  SelfAdjointEigenSolver<CXMatrix> cEigen;
  for (int i = 0; i < hull.size() - 1; i++) {
    Complex diff = hull[i + 1] - hull[i];
    Complex mul = Complex(0, 1) * diff / abs(diff);
    Value expected = mul.real() * hull[i].real() + mul.imag() * hull[i].imag();
    mulMatrix = mul * cIMat;
    testMatrix = (mulMatrix + mulMatrix.adjoint()) / 2;
    cEigen.compute(testMatrix.template selfadjointView<Lower>(), false);
    if (cEigen.info() != Success) {
      (*errStream) << line << endl;
      numErrs++;
      return false;
    }
    const auto eVals_test = cEigen.eigenvalues();
    Value edge = *max_element(eVals_test.cbegin(), eVals_test.cend());
    if (edge > expected + eps2)
      return false;
  }
  return true;
}

template <int MSize, typename Derived>
bool rnrIsPolygon(const MatrixBase<Derived>& iMat, Value eps1, Value eps2) {
  Matrix<Value, MSize, 1> cVec;
  Matrix<Value, 1, MSize - 1> rVec;
  cVec.setZero();
  rVec.setZero();
  Matrix<Value, MSize, MSize - 1> midMat;
  Matrix<Value, MSize - 1, MSize - 1> newMat;
  Value mult;
  mp_for_each<mp_iota_c<MSize - 1> > ( [&](auto I) {
    cVec += iMat.col(I);
    mult = 1 / sqrt((Value)((I + 1) * (I + 2)));
    midMat.col(I) = mult * cVec - ((I + 1) * mult) * iMat.col(I + 1);
  });
  mp_for_each<mp_iota_c<MSize - 1> > ( [&](auto I) {
    rVec += midMat.row(I);
    mult = 1 / sqrt((Value)((I + 1) * (I + 2)));
    newMat.row(I) = mult * rVec - ((I + 1) * mult) * midMat.row(I + 1);
  });
  return nrIsPolygon(newMat, eps1, eps2);
}

#define ADAPT_VALUE(z, n, unused)                                       \
  case n:                                                               \
  return rnrIsPolygon<n>(mat.template topLeftCorner<n, n>(), eps1, eps2);

template <typename Derived>
bool autoRnrIsPolygon(const MatrixBase<Derived>& mat, Value eps1, Value eps2) {
  switch (mat.cols()) {
    BOOST_PP_REPEAT_FROM_TO(MIN_VERTS, BOOST_PP_ADD(MAX_VERTS, 1), ADAPT_VALUE, 0);
  default:
    abort();
  }
}

#ifdef SELF_TEST
void lldb_hack() {
  static_vector<Complex, MAX_VERTS> a;
  static_vector<Value, MAX_VERTS> b;
  a[0] = 0;
  b[0] = 0;
}

#define TEST_MATRIX(correctAnswer, comment, ...)            \
  test << __VA_ARGS__;                                      \
  result = autoRnrIsPolygon(test, EPS1, EPS2);              \
  testCount++;                                              \
  cout << (result == correctAnswer ? "" : "not ")           \
  << "ok " << testCount << " - " << comment                 \
  << " (got " << result << "; expected "                    \
  << correctAnswer << ")" << endl

int main() {
  Matrix<Value, 6, 6> test;
  bool result;
  int testCount = 0;
  TEST_MATRIX(true, "normal (directed cycle)",
              +1,-1, 0, 0, 0, 0,
              +0, 1,-1, 0, 0, 0,
              +0, 0, 1,-1, 0, 0,
              +0, 0, 0, 1,-1, 0,
              +0, 0, 0, 0, 1,-1,
              -1, 0, 0, 0, 0, 1);
  TEST_MATRIX(true, "restricted normal",
              +4, 0,-1,-1,-1,-1,
              -1, 4, 0,-1,-1,-1,
              +0,-1, 4,-1,-1,-1,
              +0, 0, 0, 1,-1, 0,
              +0, 0, 0,-1, 1, 0,
              +0, 0, 0, 0, 0, 0);
  TEST_MATRIX(true, "pseudo normal",
              +3,-1, 0, 0,-1,-1,
              -1, 3, 0, 0,-1,-1,
              +0, 0, 1, 0,-1, 0,
              +0, 0, 0, 1, 0,-1,
              +0, 0, 0,-1, 1, 0,
              +0, 0,-1, 0, 0, 1);
  TEST_MATRIX(false, "balanced non-normal",
              +2,-1, 0,-1, 0, 0,
              +0, 1, 0, 0,-1, 0,
              -1, 0, 2, 0,-1, 0,
              +0, 0, 0, 1, 0,-1,
              -1, 0,-1, 0, 2, 0,
              +0, 0,-1, 0, 0, 1);
  cout << "1.." << testCount << endl;
  return 0;
}
#else
#include <fstream>

#ifdef NO_DYNAMIC_SIZE
typedef Matrix<Value, MAX_VERTS, MAX_VERTS> MatrixXv;
#else
typedef Matrix<Value, Dynamic, Dynamic, 0, MAX_VERTS, MAX_VERTS> MatrixXv;
#endif

MatrixXv readg_from_string(const string& iLine) {
  string line = iLine.substr(1);
  int i, j, n, c, b, sl;
#ifndef NO_DYNAMIC_SIZE
  if (line[0] == '~') {
    n = 0;
    if (line[1] == '~')
      for (c = 2; c < 8; c++)
	n = 64 * n + (line[c] - 63);
    else
      for (c = 1; c < 4; c++)
	n = 64 * n + (line[c] - 63);
  } else {
    c = 1;
    n = line[0] - 63;
  }
  MatrixXv result(n, n);
#else
  MatrixXv result;
  c = 1;
  n = MAX_VERTS;
#endif
  result.setZero();
  sl = line.size();
  for (i = c; i < sl; i++)
    line[i] -= 63;
  b = 6;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      b--;
      if (b < 0) {
	b = 5;
	c++;
      }
      result(i, j) = -(line[c] >> b) & 1;
    }
  }
  result.diagonal() = -result.rowwise().sum();
  return result;
}

int main(int argc, char** argv) {
  MatrixXv data;
  ifstream inFile;
  ofstream outFile, errFile;
  unsigned long numProcessed;
  bool useStdIn = true, useStdOut = true, useStdErr = true;
  if (argc > 1) {
    if (strlen(argv[1]) != 1 || argv[1][0] != '-') {
      inFile.open(argv[1]);
      useStdIn = false;
    }
  }
  if (argc > 2) {
    if (strlen(argv[2]) != 1 || argv[2][0] != '-') {
      outFile.open(argv[2]);
      useStdOut = false;
    }
  }
  if (argc > 3) {
    if (strlen(argv[3]) != 1 || argv[3][0] != '-') {
      errFile.open(argv[3]);
      useStdErr = false;
    }
  }
  istream &f = (useStdIn ? cin : inFile);
  ostream &o = (useStdOut ? cout : outFile);
  errStream = &(useStdErr ? cerr : errFile);
  numErrs = numProcessed = 0;
  while (getline(f, line)) {
    data = readg_from_string(line);
    numProcessed++;
#ifndef NO_DYNAMIC_SIZE
    if (autoRnrIsPolygon(data, EPS1, EPS2))
      o << line << endl;
#else
    if (rnrIsPolygon<MAX_VERTS, MatrixXv>(data, EPS1, EPS2))
      o << line << endl;
#endif
  }
  inFile.close();
  outFile.close();
  errFile.close();
  cerr << ">Z " << numProcessed << " processed, " << numErrs << " eigenvalue failures" << endl;
  return 0;
}

#endif
