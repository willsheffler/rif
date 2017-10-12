#include <Eigen/Dense>
#include <iostream>

namespace rif {
namespace cluster {

using namespace Eigen;

template <class T>
using EigV = Matrix<T, Dynamic, 1>;

template <class T>
struct Metric {
  virtual double operator()(EigV<T> a, EigV<T> b) const = 0;
  virtual ~Metric() {}
};
template <class T>
struct L2 : public Metric<T> {
  double operator()(EigV<T> a, EigV<T> b) const override {
    return (a - b).norm();
  }
};
template <class T>
struct Hamming : public Metric<T> {
  double operator()(EigV<T> a, EigV<T> b) const override {
    return (a.array() != b.array()).sum() / a.size();
  }
};
template <class T>
struct Ndiff : public Metric<T> {
  double operator()(EigV<T> a, EigV<T> b) const override {
    return (a.array() != b.array()).template cast<T>().sum();
  }
};

template <class T>
std::vector<size_t> cookie_cutter_update(
    Ref<const Matrix<T, Dynamic, Dynamic, RowMajor>> data, double thresh,
    std::vector<size_t>& centers, std::string metricstr) {
  // std::cout << "foo " << data.rows() << " " << data.cols() << std::endl;
  Metric<T>* metricp;
  if (metricstr == "L2") metricp = new L2<T>;
  if (metricstr == "hamming") metricp = new Hamming<T>;
  if (metricstr == "ndiff") metricp = new Ndiff<T>;
  Metric<T>& metric(*metricp);
  for (int i = 0; i < data.rows(); ++i) {
    bool dup = false;
    for (auto j : centers) {
      // std::cout << i << " " << j << " " << metric(data.row(i), data.row(j))
      // << std::endl;
      if (metric(data.row(i), data.row(j)) <= thresh) {
        dup = true;
        break;
      }
    }
    if (!dup) centers.push_back(i);
  }
  delete metricp;
  return centers;
  // EigV<size_t> clust(centers.size());
  // for (int i = 0; i < clust.size(); ++i) clust[i] = centers[i];
  // return clust;
}

template <class T>
std::vector<size_t> cookie_cutter(
    Ref<const Matrix<T, Dynamic, Dynamic, RowMajor>> data, double thresh,
    std::string metricstr) {
  std::vector<size_t> centers;
  return cookie_cutter_update(data, thresh, centers, metricstr);
}

//
}
}