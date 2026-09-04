#pragma once

#include "../support/checks.hpp"

#include "blocking_analysis/blocking_analysis.hpp"

#include <stdexcept>

TEST_CASE("Blocking analysis requires two complete blocks", "[blocking]") {
  BlockingAnalysis blocking{3uz};

  for (const auto energy : {1.0_fp, 2.0_fp, 3.0_fp}) {
    blocking.add(energy);
  }

  REQUIRE_FALSE(blocking.ready());
  REQUIRE_THROWS_AS(blocking.mean_and_standard_error(), std::runtime_error);

  for (const auto energy : {4.0_fp, 5.0_fp, 6.0_fp}) {
    blocking.add(energy);
  }

  REQUIRE(blocking.ready());
}

TEST_CASE("Blocking analysis uses complete block means", "[blocking]") {
  BlockingAnalysis blocking{2uz};

  for (const auto energy : {1.0_fp, 3.0_fp, 5.0_fp, 7.0_fp, 100.0_fp}) {
    blocking.add(energy);
  }

  const auto [mean, standard_error]{blocking.mean_and_standard_error()};
  require_near(mean, 4.0_fp);
  require_near(standard_error, 2.0_fp);
}
