// Copyright 2010 Martin C. Frith
#include "options.h"

#include <unistd.h>

#include <limits.h>

#include <cstdlib>  // EXIT_SUCCESS
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <getopt.h>

typedef std::runtime_error Error;

std::istream &operator>>(std::istream &s, Options::OutputType &x) {
  int i = 0;
  s >> i;
  if (i < 0 || i > 4)
    s.setstate(std::ios::failbit);
  if (s)
    x = static_cast<Options::OutputType>(i);
  return s;
}

int intFromText(const char *text) {
  long x = strtol(text, 0, 0);
  if (x > INT_MAX || x < INT_MIN) return -1;
  return x;
}

Options::Options() :
  error_rate(0.05),
  verbose(false),
  input_filename(nullptr),
  output_filename(nullptr),
  kmer_size(30),
  N_threshold(0.5),
  min_anchor_size(3),
  max_anchor_size(20),
  expand_window_size(10),
  min_copy_number(5),
  max_copy_number(1000),
  outputType(maskOut) {}

void Options::fromArgs(int argc, char **argv) {
  std::string help = "currently nothing";
  // -k for transition cost?

  const char sOpts[] = "vi:o:k:e:hn:";

  static struct option lOpts[] = {{"help", no_argument, 0, 'h'},
                                  {"verbose", no_argument, 0, 'v'},
                                  {"kmer", required_argument, 0, 'k'},
                                  {"error", required_argument, 0, 'e'},
                                  {"nth", required_argument, 0, 'n'},
                                  {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'v':
      verbose = true;
      break;
    case 'i':
      input_filename = optarg;
      break;
    case 'o':
      output_filename = optarg;
      break;
    case 'k':
      kmer_size = intFromText(optarg);
      break;
    case 'e':
      error_rate = strtod(optarg, 0);
      break;
    case 'n':
      N_threshold = intFromText(optarg);
      break;
    }
  }
}