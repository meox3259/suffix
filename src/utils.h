#pragma once

#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "zstream.h"

struct Macrosatellite {
  int left_boundary;
  int right_boundary;
  int repeat_length;
  double copy_number;
  int approximate_repeat_unit;
  double avg_match_ratio;
  std::string consensus;
  std::vector<int> anchors;

  ~Macrosatellite() {}
};

bool filter_by_seed(const std::string &s, int pos);
char safeNumberToDnaChar(int num, char default_char = 'N');
uint8_t *alloc_uint8_t(const std::string &sequence);
std::vector<std::string> split(const std::string &s, char delimiter);
bool check_substring(const char* text, const std::string& substring);

// open an input file, but if the name is "-", just return cin
std::istream &openIn(const std::string &fileName, izstream &z);

// open an output file, but if the name is "-", just return cout
std::ostream &openOut(const std::string &fileName, std::ofstream &ofs);

template <typename T> std::string stringify(const T &x) {
  std::ostringstream oss;
  oss << x;
  assert(oss);
  return oss.str();
}

template <typename T> void unstringify(T &x, const std::string &s) {
  std::istringstream iss(s);
  if (!(iss >> x) || !(iss >> std::ws).eof())
    throw std::runtime_error("can't interpret: " + s);
}

template <typename T> void unfilify(T &x, const std::string &fileName) {
  izstream z;
  std::istream &input = openIn(fileName, z);
  input >> x;
  if (!input)
    throw std::runtime_error("can't read file: " + fileName);
  // check for junk at end of file?
}

int overlap_segment(int l1, int r1, int l2, int r2);