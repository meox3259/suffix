// Copyright 2010 Martin C. Frith

#pragma once

struct Options {
  Options();

  void fromArgs(int argc, char** argv);

  bool verbose;
  char *input_filename;
  char *output_filename;
  int kmer_size;
  double error_rate;
  double N_threshold;
  int min_anchor_size;
  int max_anchor_size;
  int min_copy_number;
  int max_copy_number;
  int expand_window_size;
  int num_threads;
  int lower_bound;
  int upper_bound;
  double kmer_rate;

  enum OutputType { maskOut, probOut, countOut, bedOut, repOut } outputType;
};