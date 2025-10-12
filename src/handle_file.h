#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#define MAX_READ_LEN 1000000000
#define MAX_ID_LEN 1000
#define BLK 4096

class Read {
public:
  int len;
  uint8_t Read[MAX_READ_LEN];
  char ID[MAX_ID_LEN];
  int is_N[MAX_READ_LEN];
};

int char2int(char c);
FILE *init_file(char *inputFile);
Read *return_read(FILE *fp);