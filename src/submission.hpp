#pragma once

#include <cstddef>
#include <vector>

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid
{
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<double> grid_;

public:
  Grid(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), grid_(std::vector<double>(rows * cols)) {}

  double &operator()(std::size_t i, std::size_t j)
  {
    return this->grid_[i * this->cols_ + j];
  }
  double operator()(std::size_t i, std::size_t j) const
  {
    return this->grid_[i * this->cols_ + j];
  }

  size_t rows() const { return this->rows_; }
  size_t cols() const { return this->cols_; }
};

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid &old_grid, Grid &new_grid)
{
  size_t rows = old_grid.rows();
  size_t cols = old_grid.cols();

// Boundary on the left and right
#pragma omp parallel for
  for (int i = 0; i < rows; ++i)
  {
    new_grid(i, 0) = old_grid(i, 0);
    new_grid(i, cols - 1) = old_grid(i, cols - 1);
  }

// Bondary on the top and bottom
#pragma omp parallel for
  for (int i = 0; i < cols; ++i)
  {
    new_grid(0, i) = old_grid(0, i);
    new_grid(rows - 1, i) = old_grid(rows - 1, i);
  }

// Go through all the cells in the interior
#pragma omp parallel for
  for (int i = 1; i < rows - 1; ++i)
  {
    for (int j = 1; j < cols - 1; ++j)
    {
      new_grid(i, j) = 0.5 * old_grid(i, j) +
                       0.125 * (old_grid(i - 1, j) + old_grid(i + 1, j) +
                                old_grid(i, j - 1) + old_grid(i, j + 1));
    }
  }
}
