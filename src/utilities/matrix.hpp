#pragma once

#include <cmath>
#include <cstddef>
#include <utility> // Added for std::swap

namespace {

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
In-Place LU with partial pivoting in LU (size N*N, row-major)
pivot is length >= N storing row permutation indices
return numbers of row swaps (parity info, if you need det sign).
*/
int lower_upper_decomp(double* lowerUpper, int* pivot, std::size_t N) {
    // Track row swaps
    int swapCount{};

    for (std::size_t row = 0; row < N; ++row)
        pivot[row] = static_cast<int>(row);

    for (std::size_t col = 0; col < N; ++col) {
        // Pivot selection
        // Find row >= col maximizing |LU(row, col)|
        std::size_t pivotRow{col};
        double maxAbs{std::abs(lowerUpper[col * N + col])};

        for (std::size_t row = col + 1; row < N; ++row) {
            const double value = std::abs(lowerUpper[row * N + col]);
            if (value > maxAbs) {
                maxAbs = value;
                pivotRow = row;
            }
        }

        // max abs = 0.0 implies the pivot column is 0 & det = 0.
        if (maxAbs == 0.0)
            continue;

        if (pivotRow != col) {
            for (std::size_t col2 = 0; col2 < N; ++col2) {
                std::swap(lowerUpper[col * N + col2], lowerUpper[pivotRow * N + col2]);
            }
            std::swap(pivot[col], pivot[pivotRow]);
            ++swapCount;
        }

        // eliminate
        const double pivotValue{lowerUpper[col * N + col]};
        for (std::size_t row = col + 1; row < N; ++row) {
            lowerUpper[row * N + col] /= pivotValue; // L (i,k)
            const double multiplier{lowerUpper[row * N + col]};
            for (std::size_t col2 = col + 1; col2 < N; ++col2) {
                lowerUpper[row * N + col2] -= multiplier * lowerUpper[col * N + col2];
            }
        }
    }

    return swapCount;
}

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
solve (P^-1)LU x = b. given combined LU and pivot permutation piv.
piv encodes the row permutation applied during LU so that
we first permute b: y = P b, then solve L z = y, then U x = z.
*/
void solve_lower_upper(const double* LU, const int* pivot, const double* b, double* x,
                       std::size_t N) {
    // Apply permutation: x = Pb
    // store y in x temporarily
    for (std::size_t row = 0; row < N; ++row) {
        const std::size_t permRow{static_cast<std::size_t>(pivot[row])};
        x[row] = b[permRow];
    }

    // forward solve: ly = Pb (L has implicit on diagonal)
    for (std::size_t row = 0; row < N; ++row) {
        double sum = x[row];
        for (std::size_t col = 0; col < row; ++col) {
            sum -= LU[row * N + col] * x[col];
        }
        x[row] = sum;
    }

    // backward solve: Ux = y
    for (std::size_t rev = 0; rev < N; ++rev) {
        const std::size_t row = N - 1 - rev;
        double sum = x[row];
        for (std::size_t col = row + 1; col < N; ++col) {
            sum -= LU[row * N + col] * x[col];
        }
        x[row] = sum / LU[row * N + row];
    }
}

// Canonical representative rule for +-n deduplication:
// The canonical form is: the first nonzero component is positive
// The zero vector (0,0,0) is its own canonical representative
bool is_canonical(int n_x, int n_y, int n_z) {
    if (n_x > 0)
        return true;
    if (n_x < 0)
        return false;
    if (n_y > 0)
        return true;
    if (n_y < 0)
        return false;
    if (n_z > 0)
        return true;
    if (n_z < 0)
        return false;

    // (0,0,0)
    return true;
}

} // namespace