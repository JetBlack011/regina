// json_test.cpp: the one JSON writer (../json.h) writes what the run
// directory's writers always wrote -- its formats are frozen.

#include <string>
#include <vector>

#include "cobound/json.h"
#include "linknaming/tests/check.h"

int main() {
  // Escaping: '"' and '\' get a backslash, a newline becomes \n, and nothing
  // else changes (cascadesearch's jsonEscape, byte for byte).
  CHECK_EQ(json::escape("plain 3_1#m3_1"), std::string("plain 3_1#m3_1"), "plain text unchanged");
  CHECK_EQ(json::escape("a\"b"), std::string("a\\\"b"), "a quote");
  CHECK_EQ(json::escape("a\\b"), std::string("a\\\\b"), "a backslash");
  CHECK_EQ(json::escape("line1\nline2"), std::string("line1\\nline2"), "a newline");
  CHECK_EQ(json::escape("tab\there"), std::string("tab\there"), "a tab is left alone (as before)");
  CHECK_EQ(json::escape("PD[X[4; 1; 3; 2]]"), std::string("PD[X[4; 1; 3; 2]]"), "a PD code");
  CHECK_EQ(json::escape(""), std::string(""), "empty");
  CHECK_EQ(json::quote("say \"hi\""), std::string("\"say \\\"hi\\\"\""), "quote() wraps escape()");

  // Arrays and matrices: no spaces.
  CHECK_EQ(json::array(std::vector<int>{1, -2, 3}), std::string("[1,-2,3]"), "an array");
  CHECK_EQ(json::array(std::vector<long>{}), std::string("[]"), "an empty array");
  CHECK_EQ(json::array(std::vector<size_t>{7}), std::string("[7]"), "one element");
  CHECK_EQ(json::matrix(std::vector<std::vector<int>>{{0, 1}, {1, 0}}),
           std::string("[[0,1],[1,0]]"), "a linking matrix");
  CHECK_EQ(json::matrix(std::vector<std::vector<int>>{}), std::string("[]"), "an empty matrix");
  CHECK_EQ(json::matrix(std::vector<std::vector<int>>{{}}), std::string("[[]]"),
           "one empty row");
  return checks::finish("json_test");
}
