#include "utils.h"

#include <cstdint>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>
#include <sys/types.h>

#define overlap_threshold 0.5

enum class GENE { A, G, T, C };

const char *seeds[] = {
    "RRRRRRRYYRYY", "RRYRRYYYYYYY", "RRRRRRYYYRYY", "RRYRRRYYYYYY",
    "RRRRRYRYYRYY", "RRYRRYRYYYYY", "RRRRRYYRYRYY", "RRYRYRRYYYYY",
    "RRRRRYYYYRYY", "RRYRRRRYYYYY", "RRRRYRYYYRYY", "RRYRRRYRYYYY",
    "RRRRYYYRRYYY", "RRRYYRRRYYYY", "RRRRYYYRYRYY", "RRYRYRRRYYYY",
    "RRRYRRRYYRYY", "RRYRRYYYRYYY", "RRRYRYRYYRYY", "RRYRRYRYRYYY",
    "RRRYRYYRRYYY", "RRRYYRRYRYYY", "RRRYRYYRYRYY", "RRYRYRRYRYYY",
    "RRRYYRRYYRYY", "RRYRRYYRRYYY", "RRRYYRYRYRYY", "RRYRYRYRRYYY",
    "RRYRRRRYYRYY", "RRYRRYYYYRYY", "RRYRRYYRYRYY", "RRYRYRRYYRYY"};

const int seed_size = 32;

bool filter_by_seed(const std::string &s, int pos) {
  for (int i = 0; i < seed_size; ++i) {
    bool match = true;
    for (int j = 0; j < 12; ++j) {
      if ((seeds[i][j] == 'R' && (s[pos + j] != '0' || s[pos + j] != '2')) ||
          (seeds[i][j] == 'Y' && (s[pos + j] != '1' || s[pos + j] != '3'))) {
        match = false;
        break;
      }
    }
    if (match) {
      return true;
    }
  }
  return false;
}

char safeNumberToDnaChar(int num, char default_char) {
  if (num >= 0 && num <= 3) {
    static const char mapping[] = {'A', 'G', 'T', 'C'};
    return mapping[num];
  }
  return default_char;
}

uint8_t *alloc_uint8_t(const std::string &sequence) {
  int len = static_cast<int>(sequence.size());
  uint8_t *data = new uint8_t[len];
  for (int i = 0; i < len; ++i) {
    data[i] = static_cast<uint8_t>(sequence[i]);
  }
  return data;
}

bool N_filter(const std::string &sequence) {
  int N_count = 0;
  for (char c : sequence) {
    N_count += (c == 'N');
  }
  return N_count > 0.5 * sequence.size();
}

template <typename T>
bool N_filter(T *sequence, int len) {
  int N_count = 0;
  for (int i = 0; i < len; ++i) {
    N_count += (sequence[i] == 'N');
  }
  return N_count > 0.5 * len;
}

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, delimiter)) {
      result.push_back(item);
    }

    return result;
}

bool check_substring(const char* text, const std::string& substring) {
    if (text == nullptr || substring.empty()) {
        return false;
    }
    return strstr(text, substring.c_str()) != nullptr;
}

int overlap_segment(int l1, int r1, int l2, int r2) {
  return std::max(0, std::min(r1, r2) - std::max(l1, l2) + 1);
}

std::string firstWord(const std::string &s) {
  std::string word;
  std::istringstream iss(s);
  iss >> word;
  return word; // might be empty
}