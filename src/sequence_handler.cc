#include "sequence_handler.h"
#include "edlib_align.h"
#include "fenwick_tree.h"
#include "handle_file.h"
#include "ksw2/ksw2.h"
#include "logger.h"
#include "options.h"
#include "radix_sort.h"
#include "spoa/include/spoa/spoa.hpp"
#include "suffix_array.h"
#include "thread_pool.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <assert.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <numeric>

const int INF = std::numeric_limits<int>::max() / 2;

int match = 1, mis = 2, gap_open = 2, gap_ext = 1;
int8_t mat[25] = {1,  -2, -2, -2, 0,  -2, 1, -2, -2, 0, -2, -2, 1,
                  -2, 0,  -2, -2, -2, 1,  0, 0,  0,  0, 0,  0};

int prime_number[] = {
    -100000000, 1009, 1061, 1109, 1163, 1213, 1277, 1327, 1381,    1429,
    1481,       1531, 1579, 1627, 1693, 1741, 1789, 1847, 1901,    1949,
    1997,       2053, 2111, 2161, 2213, 2267, 2333, 2381, 2437,    2503,
    2551,       2609, 2657, 2707, 2767, 2819, 2879, 2927, 2999,    3049,
    3109,       3163, 3217, 3271, 3319, 3371, 3433, 3491, 3539,    3593,
    3643,       3691, 3739, 3793, 3847, 3907, 3967, 4000, 4050,    4100,
    4150,       4200, 4350, 4300, 4350, 4400, 4450, 4500, 4550,    4600,
    4650,       4700, 4750, 4800, 4850, 4900, 4950, 5000, +1000000};

namespace factor {

const double error_rate = 0.01;
const double kmer_rate = 0.05; // change
const double alignment_error_rate = 0.2;
const int lower_bound = 1000, upper_bound = 5000;
const int gap_delta_threshold = 50;
const int min_anchor_size = 3, max_anchor_size = 20;
const int kmer_size = 20;
const double N_threshold = 0.5;

bool is_gap_valid(int gap) { return gap >= lower_bound && gap <= upper_bound; }

} // namespace factor

namespace alignment {

std::vector<int> ksw2_get_xid(uint32_t *cigar, int n_cigar,
                              const uint8_t *query, const uint8_t *target) {
  int qi = 0;
  int ti = 0;
  std::vector<int> xid(4);
  for (int i = 0; i < n_cigar; ++i) {
    int op = cigar[i] & 0xf;
    int len = cigar[i] >> 4;
    if (op == 0) {
      // match/mismatch
      for (int j = 0; j < len; ++j) {
        if (query[qi + j] == target[ti + j]) {
          xid[0]++;
        } else {
          xid[3]++;
        }
      }
      qi += len;
      ti += len;
    } else if (op == 1) {
      xid[1] += len;
      qi += len;
    } else if (op == 2) {
      xid[2] += len;
      ti += len;
    } else {
      return {};
    }
  }
  return xid;
}

std::pair<int, int> alignment(uint8_t *query, int qlen, uint8_t *target,
                              int tlen) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  //  for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
  //    printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
  //  putchar('\n');
  fflush(stdout);
  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  return {xid[0], accumulate(xid.begin(), xid.end(), 0)};
}

std::pair<int, int> alignment(const std::string &query_string,
                              const std::string &target_string) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  uint8_t *query = alloc_uint8_t(query_string);
  uint8_t *target = alloc_uint8_t(target_string);
  int qlen = static_cast<int>(query_string.size());
  int tlen = static_cast<int>(target_string.size());
  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);

  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  delete[] query;
  delete[] target;
  return {xid[0], accumulate(xid.begin(), xid.end(), 0)};
}

int ksw2_backtrack_right_end(int n_cigar, uint32_t *cigar, uint8_t *query,
                             int qlen, uint8_t *target, int tlen,
                             double ratio) {
  int qi, ti;
  qi = ti = 0;
  std::vector<std::pair<int, int>> match;
  match.reserve(std::max(qlen, tlen));
  for (int i = 0; i < n_cigar; ++i) {
    int op = cigar[i] & 0xf;
    int len = cigar[i] >> 4;
    if (op == 0) {
      // match/mismatch
      for (int j = 0; j < len; ++j) {
        if (query[qi + j] == target[ti + j]) {
          match.push_back({1, 1});
        } else {
          match.push_back({1, 0});
        }
      }
      qi += len;
      ti += len;
    } else if (op == 1) {
      match.push_back({len, 0});
      qi += len;
    } else if (op == 2) {
      match.push_back({len, 0});
      ti += len;
    } else {
      return -1;
    }
  }

  double prefix_match_sum = 0;
  double prefix_total_sum = 0;
  int last_pos = -1;
  for (int i = 0; i < match.size(); ++i) {
    prefix_total_sum += match[i].first;
    prefix_match_sum += match[i].second;
    double curr_match_ratio = prefix_match_sum / prefix_total_sum;
    if (curr_match_ratio >= ratio) {
      last_pos = i;
    }
  }

  return last_pos;
}

int ksw2_backtrack_left_end(int n_cigar, uint32_t *cigar, uint8_t *query,
                            int qlen, uint8_t *target, int tlen, double ratio) {
  int qi, ti;
  qi = ti = 0;
  std::vector<std::pair<int, int>> match;
  match.reserve(std::max(qlen, tlen));
  for (int i = 0; i < n_cigar; ++i) {
    int op = cigar[i] & 0xf;
    int len = cigar[i] >> 4;
    if (op == 0) {
      // match/mismatch
      for (int j = 0; j < len; ++j) {
        if (query[qi + j] == target[ti + j]) {
          match.push_back({1, 1});
        } else {
          match.push_back({1, 0});
        }
      }
      qi += len;
      ti += len;
    } else if (op == 1) {
      match.push_back({len, 0});
      qi += len;
    } else if (op == 2) {
      match.push_back({len, 0});
      ti += len;
    } else {
      return -1;
    }
  }

  double prefix_match_sum = 0;
  double prefix_total_sum = 0;
  int last_pos = -1;
  for (int i = 0; i < match.size(); ++i) {
    prefix_total_sum += match[i].first;
    prefix_match_sum += match[i].second;
    double curr_match_ratio = prefix_match_sum / prefix_total_sum;
    if (curr_match_ratio >= ratio) {
      last_pos = i;
    }
  }

  return last_pos;
}

int extend_left_boundary(uint8_t *query, int qlen, uint8_t *target, int tlen,
                         double ratio) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ez));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  // 堆上分配反转缓冲
  uint8_t *rq = (uint8_t *)malloc((size_t)qlen);
  uint8_t *rt = (uint8_t *)malloc((size_t)tlen);
  if (!rq || !rt) {
    free(rq);
    free(rt);
    return -1;
  }

  for (int i = 0, j = qlen - 1; i < qlen; ++i, --j)
    rq[i] = query[j];
  for (int i = 0, j = tlen - 1; i < tlen; ++i, --j)
    rt[i] = target[j];

  ksw_extz2_sse(0, qlen, rq, tlen, rt, 5, mat, gap_open, gap_ext, w, zdrop,
                end_bonus, flag, &ez);
  int ed =
      ksw2_backtrack_left_end(ez.n_cigar, ez.cigar, rq, qlen, rt, tlen, ratio);

  free(rq);
  free(rt);
  if (ez.cigar)
    free(ez.cigar);

  return ed;
}

int extend_right_boundary(uint8_t *query, int qlen, uint8_t *target, int tlen,
                          double ratio) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  return ksw2_backtrack_right_end(ez.n_cigar, ez.cigar, query, qlen, target,
                                  tlen, ratio);
}

int ksw2_global_with_cigar(const uint8_t *query, int qlen,
                           const uint8_t *target, int tlen, int *n_cigar,
                           uint32_t **cigar) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;
  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
#ifdef __DEBUG__
  print_cigar(ez.n_cigar, ez.cigar);
#endif
  vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  int iden_n = 0;
  if (!xid.empty()) {
    iden_n = xid[0];
  }
  *n_cigar = ez.n_cigar;
  *cigar = ez.cigar;
  // if (ez.cigar) free(ez.cigar);
  return iden_n;
}

} // namespace alignment

void solve(std::ofstream &ofs, Read *Read, const Options &options) {
  LOG << "Start to build suffix array: size = " << Read->len;
  std::vector<int> sequence;
  for (int i = 0; i < Read->len; ++i) {
    sequence.push_back(Read->Read[i]);
  }

  auto suffix_array = yosupo::suffix_array(sequence);
  auto lcp_array = yosupo::lcp_array(sequence, suffix_array);

  LOG << "End to build suffix array";

  std::vector<std::vector<int>> index_set;
  std::vector<int> gap_values;
  std::vector<std::pair<int, std::vector<int>>> current_suspect_region;
  std::vector<int> right_index{suffix_array[0]};

  LOG << "Start to build suspect region";
  std::vector<int> prefix_sum_N(Read->len + 1);
  for (int i = 0; i < Read->len; ++i) {
    prefix_sum_N[i + 1] = prefix_sum_N[i] + Read->is_N[i];
  }
  //  std::cerr << "prefix_sum_N = " << prefix_sum_N[Read->len - 1] <<
  //  std::endl;
  LOG << "prefix_sum_N_last = " << prefix_sum_N[Read->len - 1];
  // 构建 suspect region

  std::vector<std::pair<int, int>> gaps;

  // 这里可以进行整体的基数排序 - M
  // 目前改进的方法为，每个k-mer匹配向后5个的所有k-mer，作为gap
  int id = 0;
  std::vector<std::pair<int, int>> right_index_pair;
  int start_position = suffix_array[0];
  int end_position = suffix_array[0];
  for (int i = 0; i < Read->len - 1; ++i) {
    if (lcp_array[i] < options.kmer_size) {
      if (right_index.size() <= factor::min_anchor_size) {
        right_index.clear();
        start_position = Read->len - 1;
        end_position = 0;
        continue;
      }
      int N_count = prefix_sum_N[end_position] - prefix_sum_N[start_position];
      if (N_count > options.N_threshold * (end_position - start_position + 1)) {
        right_index.clear();
        start_position = Read->len - 1;
        end_position = 0;
        continue;
      }
      for (int j : right_index) {
        right_index_pair.emplace_back(j, id);
      }
      id++;
      right_index.clear();
      start_position = Read->len - 1;
      end_position = 0;
    }
    right_index.push_back(suffix_array[i + 1]);
    start_position = min(start_position, suffix_array[i + 1]);
    end_position = max(end_position, suffix_array[i + 1]);
  }

  if (right_index.size() > factor::min_anchor_size) {
    int start_position = *min_element(right_index.begin(), right_index.end());
    int end_position = *max_element(right_index.begin(), right_index.end());
    int N_count = prefix_sum_N[end_position + 1] - prefix_sum_N[start_position];
    if (!(N_count >
          options.N_threshold * (end_position - start_position + 1))) {
      for (int j : right_index) {
        right_index_pair.emplace_back(j, id);
      }
      id++;
    }
  }

  LOG << "first id = " << id;
  radix_sort(
      right_index_pair, [](const auto &item) -> int { return item.first; },
      Read->len);
  LOG << "out of radix sort";

  std::vector<std::vector<int>> right_index_2d(id + 1);
  for (const auto &[index, id] : right_index_pair) {
    right_index_2d[id].push_back(index);
  }

  LOG << "right_index_2d.size() = " << right_index_2d.size();

  right_index = {suffix_array[0]};
  id = 0;
  start_position = suffix_array[0];
  end_position = suffix_array[0];
  for (int i = 0; i < Read->len - 1; ++i) {
    if (lcp_array[i] < options.kmer_size) {
      if (right_index.size() <= factor::min_anchor_size) {
        right_index.clear();
        start_position = Read->len - 1;
        end_position = 0;
        continue;
      }
      // 这里非常慢
      assert(start_position <= end_position);
      const auto &sorted_right_index = right_index_2d[id];
      sort(right_index.begin(), right_index.end());
      for (int j = 0; j < right_index.size(); ++j) {
        for (int k = 1; k <= 40 && j + k < right_index.size(); ++k) {
          int gap = right_index[j + k] - right_index[j];
          int N_count =
              prefix_sum_N[right_index[j + k]] - prefix_sum_N[right_index[j]];
          if (N_count >
              options.N_threshold * (right_index[j + k] - right_index[j] + 1)) {
            continue;
          }
          if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
            gap_values.push_back(gap);
            gaps.emplace_back(gap, right_index[j]);
          }
        }
      } /*
      for (int j = 0; j < sorted_right_index.size(); ++j) {
        for (int k = 1; k <= 40 && j + k < sorted_right_index.size(); ++k) {
          int gap = sorted_right_index[j + k] - sorted_right_index[j];
          if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
            gap_values.push_back(gap);
            gaps.emplace_back(gap, sorted_right_index[j]);
          }
        }
      } */
      id++;
      right_index.clear();
      start_position = Read->len - 1;
      end_position = 0;
    }
    right_index.push_back(suffix_array[i + 1]);
    start_position = min(start_position, suffix_array[i + 1]);
    end_position = max(end_position, suffix_array[i + 1]);
  }

  LOG << "final id = " << id;

  if (right_index.size() > factor::min_anchor_size) {
    //  auto sorted_right_index = right_index_2d[id];
    //  assert(
    //      is_sorted(sorted_right_index.begin(), sorted_right_index.end(),
    //                [](const auto &a, const auto &b) -> bool { return a < b;
    //                }));
    sort(right_index.begin(), right_index.end());
    for (int j = 0; j < right_index.size(); ++j) {
      for (int k = 1; k <= 40 && j + k < right_index.size(); ++k) {
        int gap = right_index[j + k] - right_index[j];
        int N_count =
            prefix_sum_N[right_index[j + k]] - prefix_sum_N[right_index[j]];
        if (N_count >
            options.N_threshold * (right_index[j + k] - right_index[j] + 1)) {
          continue;
        }
        if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
          gap_values.push_back(gap);
          gaps.emplace_back(gap, right_index[j]);
        }
      }
    } /*
      for (int j = 0; j < sorted_right_index.size(); ++j) {
        for (int k = 1; k <= 40 && j + k < sorted_right_index.size(); ++k) {
          int gap = sorted_right_index[j + k] - sorted_right_index[j];
          if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
            gap_values.push_back(gap);
            gaps.emplace_back(gap, sorted_right_index[j]);
          }
        }
      } */
    id++;
  }

  LOG << "gaps.size() = " << gaps.size();

  std::function<int(const std::pair<int, std::vector<int>> &)> get_index =
      [](const std::pair<int, std::vector<int>> &item) -> int {
    return item.first;
  };

  radix_sort(
      gaps, [](const auto &item) -> int { return item.second; }, Read->len);
  radix_sort(
      gap_values, [](const auto &item) -> int { return item; }, Read->len);
  gap_values.erase(std::unique(gap_values.begin(), gap_values.end()),
                   gap_values.end());

  std::vector<std::vector<pair<int, int>>> bucket_of_index(78);

  LOG << "gap_values.size() = " << gap_values.size();

  // 对gaps处理一下，挑出一些质数，然后处理
  // 这里对gap排序了，所以之后bucket_of_index内部不要排序
  for (int i = 0; i < gaps.size(); ++i) {
    int start_position = gaps[i].second;
    int end_position = start_position + gaps[i].first;
    int gap = gaps[i].first;
    int gap_index =
        std::lower_bound(prime_number, prime_number + 78, gap) - prime_number;
    int dis1 = std::abs(gap - prime_number[gap_index]);
    int dis2 = std::abs(gap - prime_number[gap_index - 1]);
    if (dis1 <= 50) {
      bucket_of_index[gap_index].emplace_back(start_position, gap);
    }
    if (dis2 <= 50) {
      bucket_of_index[gap_index - 1].emplace_back(start_position, gap);
    }
  }

  LOG << "bucket_of_index.size() = " << bucket_of_index.size();
  //  std::cerr << "bucket_of_index.size() = " << bucket_of_index.size()
  //            << std::endl;
  // 校验一下
  for (auto &vec : bucket_of_index) {
    assert(std::is_sorted(vec.begin(), vec.end(),
                          [](const auto &a, const auto &b) -> bool {
                            return a.first < b.first;
                          }));
  }

  LOG << "gap_values.size() = " << gap_values.size();
  //  std::cerr << "gap_values.size() = " << gap_values.size() << std::endl;
  std::vector<std::tuple<int, int, int>> raw_estimate_unit_region;

  // 排好序后，计算原始的
  for (int i = 0; i < bucket_of_index.size(); ++i) {
    int len_of_bucket = bucket_of_index[i].size();
    int k = 0;
    int num_of_index = 1;
    for (int j = 0; j < len_of_bucket; ++j) {
      int start_position = bucket_of_index[i][j].first;
      int end_position = start_position + bucket_of_index[i][j].second;
      int gap = bucket_of_index[i][j].second;
      num_of_index--;
      while (k < len_of_bucket && bucket_of_index[i][k].first < end_position) {
        k++;
        num_of_index++;
      }
      if (start_position > 133600000) {
        LOG << "start_position = " << start_position << ", gap = " << gap
            << " num_of_index = " << num_of_index;
      }
      if (num_of_index >= gap * factor::kmer_rate) {
        raw_estimate_unit_region.emplace_back(start_position, gap, i);
      }
    }
  }

  LOG << "raw_estimate_unit_region.size() = "
      << raw_estimate_unit_region.size();

  // 按照gap从小到大排序去做
  radix_sort(
      raw_estimate_unit_region,
      [](const auto &item) -> int { return std::get<1>(item); }, Read->len);

  LOG << "radix_sort done";

  // 树状数组
  OptimizedMaxFenwickTree ft(Read->len + 1);

  int seq_len = Read->len;
  auto already_covered = [&ft, &seq_len](int start_position, int end_position,
                                         int unit_size) {
    int l_max = ft.queryPrefixMax(min(start_position, seq_len));
    return l_max >= end_position;
  };

  char *query = new char[seq_len];
  for (int i = 0; i < seq_len; ++i) {
    query[i] = safeNumberToDnaChar(Read->Read[i]);
  }

  int cnt = 0;
  std::vector<Macrosatellite> macrosatellites;
  // 现在是按照unit_size排序
  // 找到当前gap对应的bucket，假设上次区间分割点为[s, e]，
  // 下次查找两个匹配的k-mer (s1, e1)，其中(s < s1 < e <
  // e1)，且s1离s和e1离e都尽量近 再做alignment，向右延伸

  ThreadPool pool(options.num_threads);
  std::vector<std::future<void>> results;
  std::mutex macrosatellites_mutex;
  for (const auto &[start_position, unit_size, index] :
       raw_estimate_unit_region) {
    if (start_position > 133600000) {
      LOG << "start_position = " << start_position
          << ", unit_size = " << unit_size;
    }
    int end_position = start_position + unit_size;
    if (already_covered(start_position + 1, end_position, unit_size)) {
      continue;
    }
    int s, e;
    s = start_position;
    e = end_position;

    auto check_valid = [&](int new_len, vector<int> &anchors) {
      sort(anchors.begin(), anchors.end());
      anchors.erase(unique(anchors.begin(), anchors.end()), anchors.end());
      for (int i = 1; i < anchors.size(); ++i) {
        int len = anchors[i] - anchors[i - 1] + 1;
        if (abs(len - new_len) > min(len, new_len) * 0.05) {
          LOG << "len = " << len << ", new_len = " << new_len;
          return false;
        }
      }
      return true;
    };

    std::vector<int> anchors{start_position, end_position};
    while (e < Read->len) {
      int qlen = unit_size;
      if (Read->len - e < qlen * 0.8) {
        break;
      }
      int tlen = unit_size * (1 + options.repeat_unit_error_rate);
      int delta = unit_size * options.repeat_unit_error_rate;
      int s1, e1;
      int l, r, l0, r0;

      l = s;
      r = e;
      l0 = e;
      r0 = e + (e - s + 1) + delta;
      qlen = r - l + 1;
      tlen = r0 - l0 + 1;

      int overlap_region = overlap_segment(l, r, l0, r0);
      if (overlap_region > unit_size * 0.8) {
        break;
      }

      int ed = edlib_align_HW(query + l, qlen, query + l0, tlen, &s1, &e1,
                              std::min(qlen, tlen));
      LOG << "ed = " << ed;
      if (ed > min(qlen, tlen) * options.extend_error_rate) {
        LOG << "ed = " << ed;
        break;
      }

      int new_len = e1 - s1 + 1;
      //  if (!check_valid(new_len, anchors)) {
      //    LOG << "111";
      //    break;
      //  }
      if (abs(new_len - unit_size) > unit_size * 0.05) {
        LOG << "111" << " " << "new_len = " << new_len << " "
            << "unit_size = " << unit_size << " ed = " << ed;
        break;
      }

      anchors.push_back(e);
      s = l0 + s1;
      e = l0 + e1;
    }

    anchors.push_back(e);

    s = start_position;
    e = end_position;

    while (s >= 0) {
      int qlen = unit_size;
      if (s < qlen * 0.8) {
        break;
      }
      int tlen = unit_size * (1 + options.repeat_unit_error_rate);
      int delta = unit_size * options.repeat_unit_error_rate;
      int s1, e1;
      int l, r, l0, r0;
      l = s;
      r = e;
      l0 = std::max(0, s - delta - tlen);
      r0 = l0 + tlen;
      qlen = r - l + 1;
      tlen = r0 - l0 + 1;

      int overlap_region = overlap_segment(l, r, l0, r0);
      if (overlap_region > unit_size * 0.8) {
        break;
      }

      int ed = edlib_align_HW(query + l, qlen, query + l0, tlen, &s1, &e1,
                              std::min(qlen, tlen));
      LOG << "ed = " << ed
          << " tolr = " << min(qlen, tlen) * options.extend_error_rate
          << "tlen = " << tlen << " qlen = " << qlen;
      if (ed > min(qlen, tlen) * options.extend_error_rate) {
        LOG << "ed = " << ed;
        break;
      }

      int new_len = e1 - s1 + 1;
      //  if (!check_valid(new_len, anchors)) {
      //    LOG << "111";
      //    break;
      //  }
      if (abs(new_len - unit_size) > unit_size * 0.05) {
        LOG << "111" << " " << "new_len = " << new_len << " "
            << "unit_size = " << unit_size << " ed = " << ed;
        break;
      }

      anchors.push_back(s);
      s = l0 + s1;
      e = l0 + e1;
    }

    anchors.push_back(s);

    sort(anchors.begin(), anchors.end());
    anchors.erase(unique(anchors.begin(), anchors.end()), anchors.end());

    int raw_left_boundary = anchors[0];
    int raw_right_boundary = anchors.back();
    int repeat_length = raw_right_boundary - raw_left_boundary + 1;

    int copy_number = static_cast<int>(anchors.size()) - 1;
    int approximate_repeat_unit = repeat_length / copy_number;
    int delta = approximate_repeat_unit / 2;

    int left_extend_end, right_extend_end;
    left_extend_end = std::max(0, anchors[0] - delta);
    right_extend_end = std::min(Read->len, anchors.back() + delta);

    LOG << "anchors.size() = " << anchors.size() << " unit_size = " << unit_size
        << " start_position = " << start_position;
    if (anchors.size() <= 5) {
      ft.updateMax(left_extend_end + 1, right_extend_end);
      continue;
    }

    int anchor_size = static_cast<int>(anchors.size());
    std::vector<std::string> sequences;
    sequences.reserve(anchor_size);
    for (int i = 0; i < anchor_size - 1; ++i) {
      std::string cur_seq = "";
      for (int j = anchors[i]; j < anchors[i + 1]; ++j) {
        cur_seq += safeNumberToDnaChar(Read->Read[j]);
      }
      sequences.push_back(cur_seq);
    }

    int l0 = std::max(0, raw_left_boundary - unit_size + 1);
    int r0 = raw_left_boundary;
    int len0 = r0 - l0 + 1;

    uint8_t *sequence_middle =
        const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(
            sequences[sequences.size() / 2].c_str()));
    int sequence_middle_len =
        static_cast<int>(sequences[sequences.size() / 2].size());

    int move_left_boundary_len =
        alignment::extend_left_boundary(sequence_middle, sequence_middle_len,
                                        reinterpret_cast<uint8_t *>(query) + l0,
                                        len0, 1. - options.error_rate * 2);

    int l1 = raw_right_boundary;
    int r1 = std::min(Read->len, raw_right_boundary + unit_size - 1);
    int len1 = r1 - l1 + 1;

    int move_right_boundary_len = alignment::extend_right_boundary(
        sequence_middle, sequence_middle_len,
        reinterpret_cast<uint8_t *>(query) + l1, len1,
        1. - options.error_rate * 2);

    LOG << "move_left_boundary_len = " << move_left_boundary_len
        << " move_right_boundary_len = " << move_right_boundary_len;
    raw_left_boundary -= move_left_boundary_len + 1;
    raw_right_boundary += move_right_boundary_len + 1;
    anchors.push_back(raw_left_boundary);
    anchors.push_back(raw_right_boundary);
    sort(anchors.begin(), anchors.end());
    anchors.erase(unique(anchors.begin(), anchors.end()), anchors.end());

    ft.updateMax(
        std::max(raw_left_boundary + 1 - approximate_repeat_unit / 2, 1),
        std::min(raw_right_boundary + approximate_repeat_unit / 2, Read->len));

    sequences.clear();
    anchor_size = static_cast<int>(anchors.size());
    for (int i = 0; i < anchor_size - 1; ++i) {
      std::string cur_seq = "";
      for (int j = anchors[i]; j < anchors[i + 1]; ++j) {
        cur_seq += safeNumberToDnaChar(Read->Read[j]);
      }
      sequences.push_back(cur_seq);
    }

    assert(sequences.size() != 0);

    results.emplace_back(pool.enqueue([&macrosatellites, &macrosatellites_mutex,
                                       sequences, anchors, raw_left_boundary,
                                       raw_right_boundary]() mutable {
      auto alignment_engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps
      spoa::Graph graph{};

      for (auto &sequence : sequences) {
        auto alignment = alignment_engine->Align(sequence, graph);
        graph.AddAlignment(alignment, sequence);
      }

      auto consensus = graph.GenerateConsensus();
      int consensus_len = static_cast<int>(consensus.size());
      uint8_t *cons =
          reinterpret_cast<uint8_t *>(const_cast<char *>(consensus.c_str()));

      Macrosatellite macrosatellite;
      macrosatellite.left_boundary = raw_left_boundary;
      macrosatellite.right_boundary = raw_right_boundary;
      macrosatellite.repeat_length = raw_right_boundary - raw_left_boundary + 1;
      macrosatellite.copy_number =
          static_cast<double>(raw_right_boundary - raw_left_boundary + 1) /
          static_cast<double>(consensus_len);
      macrosatellite.approximate_repeat_unit = consensus_len;
      macrosatellite.consensus = std::move(consensus);
      macrosatellite.anchors = std::move(anchors);

      double total_match = 0.;
      double total_length = 0.;

      bool ok = true;
      for (int i = 1; i + 1 < sequences.size(); ++i) {
        int sz = sequences[i].size();
        int delta = abs(consensus_len - sz);
        if (delta >= consensus_len * 0.2) {
          ok = false;
        }
      }
      if (!ok) {
        LOG << "out of range";
        return;
      }
      int cnt = 0;
      double avg_match_ratio = 0;
      for (auto &sequence : sequences) {
        const uint8_t *seq =
            reinterpret_cast<const uint8_t *>(sequence.c_str());
        int seq_len = static_cast<int>(sequence.size());
        auto [match, total] = alignment::alignment(
            const_cast<uint8_t *>(seq), seq_len, cons, consensus_len);
        cnt++;
        avg_match_ratio +=
            static_cast<double>(match) / static_cast<double>(sequence.size());
        LOG << "match = " << match << " sizes = " << sequence.size()
            << " total = " << total;
      }
      macrosatellite.avg_match_ratio = avg_match_ratio / cnt;
      {
        std::lock_guard<std::mutex> lock(macrosatellites_mutex);
        macrosatellites.emplace_back(macrosatellite);
      }
    }));
  }

  for (auto &result : results) {
    result.get();
  }

  LOG << "cnt = " << cnt;
  LOG << "macrosatellites.size() = " << macrosatellites.size();

  std::string name = firstWord(Read->ID);

  sort(macrosatellites.begin(), macrosatellites.end(),
       [](const Macrosatellite &lhs, const Macrosatellite &rhs) -> bool {
         return lhs.left_boundary < rhs.left_boundary;
       });

  std::vector<bool> mark(macrosatellites.size(), false);
  for (int i = 0; i < macrosatellites.size(); ++i) {
    for (int j = 0; j < macrosatellites.size(); ++j) {
      if (i == j)
        continue;
      if (macrosatellites[i].left_boundary >=
              macrosatellites[j].left_boundary &&
          macrosatellites[i].right_boundary <=
              macrosatellites[j].right_boundary) {
        mark[i] = true;
        break;
      }
    }
  }
  for (int i = 0; i < macrosatellites.size(); ++i) {
    auto macrosatellite = macrosatellites[i];
    if (macrosatellite.copy_number <= 5)
      continue;
    if (mark[i])
      continue;
    std::cout << name << ' ' << macrosatellite.left_boundary << ' '
              << macrosatellite.right_boundary << ' '
              << macrosatellite.approximate_repeat_unit << ' '
              << macrosatellite.copy_number << ' '
              << macrosatellite.avg_match_ratio << ' '
              << macrosatellite.consensus << ' ';

    for (int i = 0; i < macrosatellite.anchors.size(); ++i) {
      std::cout << macrosatellite.anchors[i]
                << ",\n"[i + 1 == macrosatellite.anchors.size()];
    }
  }
}