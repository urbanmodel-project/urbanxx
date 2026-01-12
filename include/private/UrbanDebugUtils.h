#ifndef URBAN_DEBUG_UTILS_H
#define URBAN_DEBUG_UTILS_H

#include <Kokkos_Core.hpp>
#include <iomanip>
#include <iostream>
#include <string>

namespace URBANXX {

// Debug utility to print 1D Kokkos::View to console
// Supports limiting output for large views
template <typename ViewType>
void print_view_1d(const ViewType &view, const std::string &name = "",
                   int precision = 15, std::size_t max_elements = 100) {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " = [";
  const std::size_t n = std::min(h_view.extent(0), max_elements);
  for (std::size_t i = 0; i < n; ++i) {
    std::cout << std::scientific << std::setprecision(precision) << h_view(i);
    if (i + 1 < n)
      std::cout << ", ";
  }
  if (h_view.extent(0) > max_elements)
    std::cout << ", ... (" << h_view.extent(0) - max_elements << " more)";
  std::cout << "]\n";
}

// Debug utility to print 2D Kokkos::View to console
// Supports limiting output for large views
template <typename ViewType>
void print_view_2d(const ViewType &view, const std::string &name = "",
                   int precision = 15, std::size_t max_rows = 50,
                   std::size_t max_cols = 50) {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " (" << h_view.extent(0) << "x" << h_view.extent(1)
            << "):\n";
  const std::size_t nrows = std::min(h_view.extent(0), max_rows);
  const std::size_t ncols = std::min(h_view.extent(1), max_cols);
  for (std::size_t i = 0; i < nrows; ++i) {
    for (std::size_t j = 0; j < ncols; ++j)
      std::cout << std::scientific << std::setprecision(precision)
                << h_view(i, j) << " ";
    if (h_view.extent(1) > max_cols)
      std::cout << "... (" << h_view.extent(1) - max_cols << " more)";
    std::cout << "\n";
  }
  if (h_view.extent(0) > max_rows)
    std::cout << "... (" << h_view.extent(0) - max_rows << " more rows)\n";
}

// Debug utility to print 3D Kokkos::View to console
// Supports limiting output for large views
template <typename ViewType>
void print_view_3d(const ViewType &view, const std::string &name = "",
                   int precision = 15, std::size_t max_i = 10,
                   std::size_t max_j = 10, std::size_t max_k = 10) {
  auto h_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  std::cout << name << " (" << h_view.extent(0) << "x" << h_view.extent(1)
            << "x" << h_view.extent(2) << "):\n";
  const std::size_t ni = std::min(h_view.extent(0), max_i);
  const std::size_t nj = std::min(h_view.extent(1), max_j);
  const std::size_t nk = std::min(h_view.extent(2), max_k);
  for (std::size_t i = 0; i < ni; ++i) {
    std::cout << "i=" << i << ":\n";
    for (std::size_t j = 0; j < nj; ++j) {
      for (std::size_t k = 0; k < nk; ++k)
        std::cout << std::scientific << std::setprecision(precision)
                  << h_view(i, j, k) << " ";
      if (h_view.extent(2) > max_k)
        std::cout << "... (" << h_view.extent(2) - max_k << " more)";
      std::cout << "\n";
    }
    if (h_view.extent(1) > max_j)
      std::cout << "... (" << h_view.extent(1) - max_j << " more rows)\n";
    std::cout << "\n";
  }
  if (h_view.extent(0) > max_i)
    std::cout << "... (" << h_view.extent(0) - max_i << " more slices)\n";
}

} // namespace URBANXX

#endif // URBAN_DEBUG_UTILS_H
