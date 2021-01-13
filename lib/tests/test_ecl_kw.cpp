
#include <catch2/catch.hpp>
#include <ert/ecl/ecl_kw.hpp>

using namespace Catch;
using namespace Matchers;

TEST_CASE("Test header initilization", "[unittest]") {
    GIVEN("A ecl_kw is created with a header longer than 8 characters"){
        std::string header ="a_keyword_header_longer_than_8";
        ecl_kw_type * kw = ecl_kw_alloc_new(header.c_str(), 3, ECL_FLOAT, NULL);
        THEN("get_header returns the same name"){
            REQUIRE(std::string(ecl_kw_get_header(kw)) == header);
        }
        THEN("get_header8 returns the first 8 characters of the name") {
            REQUIRE(std::string(ecl_kw_get_header8(kw)) == header.substr(0,8));
        }

        WHEN("The header name is set"){
            std::string other_header ="some_other_name";
            ecl_kw_set_header_name(kw, other_header.c_str());

            THEN("get_header returns the same name"){
                REQUIRE(std::string(ecl_kw_get_header(kw)) == other_header);
            }
            THEN("get_header8 returns the first 8 characters of the name") {
                REQUIRE(std::string(ecl_kw_get_header8(kw)) == other_header.substr(0,8));
            }
        }
        ecl_kw_free(kw);
    }
    GIVEN("A ecl_kw is created with a header less than 8 characters"){
        std::string header ="PORO";
        ecl_kw_type * kw = ecl_kw_alloc_new(header.c_str(), 3, ECL_FLOAT, NULL);
        THEN("get_header returns the same name"){
            REQUIRE(std::string(ecl_kw_get_header(kw)) == header);
        }
        THEN("get_header8 returns a padded name") {
            REQUIRE(std::string(ecl_kw_get_header8(kw)) == "PORO    ");
        }
        ecl_kw_free(kw);
    }
}
