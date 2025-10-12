#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "handle_file.h"
#include "logger.h"
#include "options.h"
#include "sequence_handler.h"
#include "utils.h"
//#include "fasta_sequence.h"

bool LogStream::verbose = false;

class Timer {
  std::chrono::steady_clock::time_point start;

public:
  Timer() : start(std::chrono::steady_clock::now()) {}
  ~Timer() {
    auto end = std::chrono::steady_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "程序运行时间: " << duration.count() << " 毫秒" << std::endl;
  }
};

int main(int argc, char *argv[]) {
  Options options;
  options.fromArgs(argc, argv);

  if (options.verbose) {
    LogStream::verbose = true;
  }

  assert(options.input_filename != nullptr);
  Timer timer;

  std::ofstream ofs(options.output_filename);

  FILE *fp = init_file(options.input_filename);
  while (true) {
    LOG << "start while";
    Read *Read = return_read(fp);
    LOG << "end return_read";
    if (Read == nullptr) {
      break;
    }
    if (Read->len == 0) {
      break;
    }
    LOG << "len = " << Read->len;
    LOG << "Read->ID = " << Read->ID;
    solve(ofs, Read, options);
    delete Read;
    LOG << "end solve";
  }
}