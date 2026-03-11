#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace slater_plane_wave_detail {

struct NVector {
    int x;
    int y;
    int z;
    int mag_sq;
};

struct ShellFillingResult {
    std::vector<NVector> unique_n_vectors;
    std::vector<std::size_t> orbital_k_index;
    std::vector<std::uint8_t> orbital_type;
};

[[nodiscard]] inline bool is_canonical(int n_x, int n_y, int n_z) noexcept {
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

    return true;
}

[[nodiscard]] inline int initial_n_max_estimate(std::size_t num_orbitals) noexcept {
    return static_cast<int>(std::ceil(std::cbrt(static_cast<double>(num_orbitals)))) + 2;
}

[[nodiscard]] inline bool radius_covers_shell(int n_max, int max_mag_sq) noexcept {
    const int max_component_needed{static_cast<int>(std::floor(std::sqrt(static_cast<double>(max_mag_sq))))};
    return n_max >= max_component_needed;
}

[[nodiscard]] inline ShellFillingResult generate_shell_filling(std::size_t num_orbitals, int initial_n_max = -1) {
    ShellFillingResult result{};
    result.orbital_k_index.resize(num_orbitals);
    result.orbital_type.resize(num_orbitals, 0);

    int n_max{initial_n_max >= 0 ? initial_n_max : initial_n_max_estimate(num_orbitals)};
    if (n_max < 0) {
        n_max = 0;
    }

    while (true) {
        const std::size_t side{static_cast<std::size_t>((2 * n_max) + 1)};

        std::vector<NVector> n_candidates{};
        n_candidates.reserve(side * side * side);

        for (int new_x = -n_max; new_x <= n_max; ++new_x) {
            for (int new_y = -n_max; new_y <= n_max; ++new_y) {
                for (int new_z = -n_max; new_z <= n_max; ++new_z) {
                    if (!is_canonical(new_x, new_y, new_z)) {
                        continue;
                    }

                    const int mag_sq{new_x * new_x + new_y * new_y + new_z * new_z};
                    n_candidates.emplace_back(new_x, new_y, new_z, mag_sq);
                }
            }
        }

        std::sort(n_candidates.begin(), n_candidates.end(), [](const NVector& a, const NVector& b) {
            if (a.mag_sq != b.mag_sq) {
                return a.mag_sq < b.mag_sq;
            }
            if (a.x != b.x) {
                return a.x < b.x;
            }
            if (a.y != b.y) {
                return a.y < b.y;
            }
            return a.z < b.z;
        });

        result.unique_n_vectors.clear();
        result.unique_n_vectors.reserve(num_orbitals);

        std::size_t orb_idx{};
        int last_mag_sq{};

        for (const NVector& candidate : n_candidates) {
            if (orb_idx >= num_orbitals) {
                break;
            }

            const std::size_t k_idx{result.unique_n_vectors.size()};
            result.unique_n_vectors.push_back(candidate);

            result.orbital_k_index[orb_idx] = k_idx;
            result.orbital_type[orb_idx] = 0;
            ++orb_idx;

            if (candidate.mag_sq != 0 && orb_idx < num_orbitals) {
                result.orbital_k_index[orb_idx] = k_idx;
                result.orbital_type[orb_idx] = 1;
                ++orb_idx;
            }

            last_mag_sq = candidate.mag_sq;
        }

        if (orb_idx == num_orbitals && radius_covers_shell(n_max, last_mag_sq)) {
            return result;
        }

        n_max = std::max(1, n_max * 2);
    }
}

} // namespace slater_plane_wave_detail
