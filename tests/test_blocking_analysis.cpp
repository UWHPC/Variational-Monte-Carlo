#include "test_utilities.hpp"

#include "blocking_analysis/blocking_analysis.hpp"

#include <stdexcept>

TEST_CASE("BlockingAnalysis readiness requires two complete blocks", "[blocking]") {
  BlockingAnalysis blocking{3U};

  REQUIRE_FALSE(blocking.ready());
  REQUIRE_THROWS_AS(blocking.mean_and_standard_error(), std::runtime_error);

  blocking.add(1.0_fp);
  blocking.add(2.0_fp);
  blocking.add(3.0_fp);

  REQUIRE_FALSE(blocking.ready());
  REQUIRE_THROWS_AS(blocking.mean_and_standard_error(), std::runtime_error);

  blocking.add(4.0_fp);
  blocking.add(5.0_fp);
  blocking.add(6.0_fp);

  REQUIRE(blocking.ready());
}

TEST_CASE("BlockingAnalysis computes mean and standard error from block means", "[blocking]") {
  BlockingAnalysis blocking{2U};
  blocking.add(1.0_fp);
  blocking.add(3.0_fp);
  blocking.add(5.0_fp);
  blocking.add(7.0_fp);

  const auto [mean, standardError]{blocking.mean_and_standard_error()};
  require_near(mean, 4.0_fp);
  require_near(standardError, 2.0_fp);
}

TEST_CASE("BlockingAnalysis ignores incomplete trailing blocks", "[blocking]") {
  BlockingAnalysis blocking{2U};
  blocking.add(1.0_fp);
  blocking.add(2.0_fp);
  blocking.add(3.0_fp);
  blocking.add(4.0_fp);
  blocking.add(100.0_fp);

  REQUIRE(blocking.ready());

  const auto [mean, standardError]{blocking.mean_and_standard_error()};
  require_near(mean, 2.5_fp);
  require_near(standardError, 1.0_fp);
}