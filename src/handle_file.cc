#include "handle_file.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

int char2int(char c) {
  if (c == 'A' || c == 'a') {
    return 0;
  }
  if (c == 'G' || c == 'g') {
    return 1;
  }
  if (c == 'T' || c == 't') {
    return 2;
  }
  if (c == 'C' || c == 'c') {
    return 3;
  }
  return -1;
}

FILE *init_file(char *inputFile) {
  FILE *fp = fopen(inputFile, "r");
  if (fp == NULL) {
    std::cerr << "fatal error: cannot open " << inputFile << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return fp;
}

Read *return_read(FILE *fp) {
  char *s = new char[BLK];
  static char *ReadID = new char[BLK];
  static int ReadsNum = 0;
  int i;
  int cnt = 0;
  bool hasRead = false;
  Read *CurrentRead = new Read;
  while (fgets(s, BLK, fp) != NULL) {
    hasRead = true;
    if (s[0] == '>') {
      if (ReadsNum == 0) {
        for (i = 1; s[i] != '\0' && s[i] != '\n' && s[i] != '\r'; ++i) {
          ReadID[i - 1] = s[i];
        }
        ReadID[i - 1] = '\0';
        ReadsNum++;
      } else {
        int j;
        for (j = 0; ReadID[j] != '\0' && ReadID[j] != '\n' && ReadID[j] != '\r';
             ++j) {
          CurrentRead->ID[j] = ReadID[j];
        }
        CurrentRead->ID[j] = '\0';
        for (i = 1; s[i] != '\0' && s[i] != '\n' && s[i] != '\r'; ++i) {
          ReadID[i - 1] = s[i];
        }
        ReadID[i - 1] = '\0';
        CurrentRead->len = cnt;
        delete[] s;
        return CurrentRead;
      }
    } else {
      for (i = 0; s[i] != '\0' && s[i] != '\n' && s[i] != '\r'; ++i) {
        if (!char2int(s[i]) == -1) {
          continue;
        }

        CurrentRead->Read[cnt++] = char2int(s[i]);
        CurrentRead->is_N[cnt - 1] = s[i] == 'N';
      }
    }
  }
  if (!hasRead) {
    CurrentRead->len = 0;
    delete[] s;
    //  delete[] ReadID;
    return CurrentRead;
  }
  int j;
  for (j = 0; ReadID[j] != '\0' && ReadID[j] != '\n' && ReadID[j] != '\r';
       ++j) {
    CurrentRead->ID[j] = ReadID[j];
  }
  CurrentRead->ID[j] = '\0';
  CurrentRead->len = cnt;
  delete[] s;
  delete[] ReadID;
  return CurrentRead;
}