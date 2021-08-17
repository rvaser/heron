// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <unordered_map>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "edlib.h"  // NOLINT
#include "ksw2.h"  // NOLINT
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"coverage", required_argument, nullptr, 'c'},
  {"low", required_argument, nullptr, 'L'},
  {"high", required_argument, nullptr, 'H'},
  {"frequencies", no_argument, nullptr, 'f'},
  {"ksw2", no_argument, nullptr, 'k'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[heron::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: heron [options ...] <sequences>\n"
      "\n"
      "  # default output is to stdout\n"
      "  <sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    --coverage <int>\n"
      "      default: 3\n"
      "      declare solid base with coverage greater than provided value\n"
      "    --low <double>\n"
      "      default: 0.333\n"
      "      declare solid base with frequency greater than provided value\n"
      "    --high <double>\n"
      "      default: 0.666\n"
      "      declare solid base with frequency less than provided value\n"
      "    -f, --frequencies\n"
      "      default: \n"
      "      technique to detect variants\n"
      "    -k, --ksw2\n"
      "      default: edlib\n"
      "      used alignment method\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::uint32_t coverage = 3;
  double low = 0.333;
  double high = 0.666;

  bool use_frequencies = false;
  bool use_ksw2 = false;
  std::uint32_t num_threads = 1;

  std::string optstr = "fkt:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'c': coverage = std::atoi(optarg); break;
      case 'L': low = std::atof(optarg); break;
      case 'H': high = std::atof(optarg); break;
      case 'f': use_frequencies = true; break;
      case 'k': use_ksw2 = true; break;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc) {
    std::cerr << "[heron::] error: missing input file!" << std::endl;
    return 1;
  }

  auto sparser = CreateParser(argv[optind]);
  if (sparser == nullptr) {
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  try {
    sequences = sparser->Parse(-1);
  } catch (const std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  std::cerr << "[heron::] parsed " << sequences.size() << " sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  struct Pile {
    std::uint32_t a;
    std::uint32_t c;
    std::uint32_t g;
    std::uint32_t t;
    std::uint32_t i;
  };

  std::vector<std::vector<Pile>> piles(sequences.size());
  for (const auto& it : sequences) {
    piles[it->id].resize(it->inflated_len);
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
  ram::MinimizerEngine minimizer_engine{thread_pool, 15U, 5U};

  auto edlib_wrapper = [&] (
      std::uint32_t i,
      const biosoup::Overlap& it,
      const std::string& lhs,
      const std::string& rhs) -> void {
    EdlibAlignResult result = edlibAlign(
        lhs.c_str(), lhs.size(),
        rhs.c_str(), rhs.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
    if (result.status == EDLIB_STATUS_OK) {
      std::uint32_t lhs_pos = it.lhs_begin;
      std::uint32_t rhs_pos = 0;
      for (int j = 0; j < result.alignmentLength; ++j) {
        switch (result.alignment[j]) {
          case 0:
          case 3: {
            switch (rhs[rhs_pos]) {
              case 'A': ++piles[i][lhs_pos].a; break;
              case 'C': ++piles[i][lhs_pos].c; break;
              case 'G': ++piles[i][lhs_pos].g; break;
              case 'T': ++piles[i][lhs_pos].t; break;
              default: break;
            }
            ++lhs_pos;
            ++rhs_pos;
            break;
          }
          case 1: {
            ++piles[i][lhs_pos].i;
            ++lhs_pos;
            break;
          }
          case 2: {
            ++rhs_pos;
            break;
          }
          default: break;
        }
      }
    }
    edlibFreeAlignResult(result);
  };

  auto ksw2_wrapper = [&] (
      std::uint32_t i,
      const biosoup::Overlap& it,
      const std::string& lhs,
      const std::string& rhs) -> void {
    std::int8_t m = 3;
    std::int8_t n = -5;
    std::int8_t g = 4;
    std::int8_t e = 4;
    std::int8_t mn[25] = {
        m, n, n, n, 0,
        n, m, n, n, 0,
        n, n, m, n, 0,
        n, n, n, m, 0,
        0, 0, 0, 0, 0
    };

    std::unordered_map<char, std::uint8_t> transform = {
      {'A', 0}, {'a', 0},
      {'C', 1}, {'c', 1},
      {'G', 2}, {'g', 2},
      {'T', 3}, {'t', 3}
    };

    auto lhs_ = new std::uint8_t[lhs.size()];
    for (std::size_t j = 0; j < lhs.size(); ++j) {
      lhs_[j] = transform[lhs[j]];
    }

    auto rhs_ = new std::uint8_t[rhs.size()];
    for (std::size_t j = 0; j < rhs.size(); ++j) {
      rhs_[j] = transform[rhs[j]];
    }

    int m_cigar = 0, n_cigar = 0;
    std::uint32_t* cigar = nullptr;

    auto score = ksw_gg2_sse(
        nullptr,   // void *km
        lhs.size(),  // int qlen
        lhs_,  // const uint8_t *query
        rhs.size(),  // int tlen
        rhs_,  // const uint8_t *target
        5,  // int8_t m
        mn,  // const int8_t *mat
        g,  // int8_t gapo
        e,  // int8_t gape
        500,  // int w
        &m_cigar,  // int *m_cigar_
        &n_cigar,  // int *n_cigar_
        &cigar);  // uint32_t **cigar_

    if (score > 0 && n_cigar > 0) {
      std::uint32_t lhs_pos = it.lhs_begin;
      std::uint32_t rhs_pos = 0;
      for (std::size_t j = 0; j < static_cast<std::size_t>(n_cigar); ++j) {
        std::size_t count = cigar[j] >> 4;
        std::size_t op = cigar[j] & 15;
        switch (op) {
          case 0: {  // M
            for (std::size_t k = 0; k < count; ++k) {
              switch (rhs[rhs_pos + k]) {
                case 'A': ++piles[i][lhs_pos + k].a; break;
                case 'C': ++piles[i][lhs_pos + k].c; break;
                case 'G': ++piles[i][lhs_pos + k].g; break;
                case 'T': ++piles[i][lhs_pos + k].t; break;
                default: break;
              }
            }
            lhs_pos += count;
            rhs_pos += count;
            break;
          }
          case 1: {  // I
            for (std::size_t k = 0; k < count; ++k) {
              ++piles[i][lhs_pos + k].i;
            }
            lhs_pos += count;
            break;
          }
          case 2: {  // D
            rhs_pos += count;
            break;
          }
          default: break;
        }
      }
    }

    free(cigar);

    delete[] rhs_;
    delete[] lhs_;
  };

  std::size_t bytes = 0;
  for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
    bytes += sequences[i]->inflated_len;
    if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
      continue;
    }
    bytes = 0;

    timer.Start();

    minimizer_engine.Minimize(
        sequences.begin() + j,
        sequences.begin() + i + 1,
        true);
    minimizer_engine.Filter(0.001);

    std::cerr << "[heron::] minimized "
              << j << " - " << i + 1 << " / " << sequences.size() << " "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    std::vector<std::future<void>> futures;
    for (std::uint32_t k = 0; k < i + 1; ++k) {
      futures.emplace_back(thread_pool->Submit(
          [&] (std::uint32_t i) -> void {
            auto overlaps = minimizer_engine.Map(sequences[i], true, false, true);  // NOLINT

            for (const auto& it : overlaps) {
              auto lhs = sequences[i]->InflateData(
                  it.lhs_begin,
                  it.lhs_end - it.lhs_begin);
              biosoup::NucleicAcid rhs_{"", sequences[it.rhs_id]->InflateData(
                  it.rhs_begin,
                  it.rhs_end - it.rhs_begin)};
              if (!it.strand) {
                rhs_.ReverseAndComplement();
              }
              auto rhs = rhs_.InflateData();

              if (use_ksw2) {
                ksw2_wrapper(i, it, lhs, rhs);
              } else {
                edlib_wrapper(i, it, lhs, rhs);
              }
            }
          },
          k));

      bytes += sequences[k]->inflated_len;
      if (k != i && bytes < (1U << 30)) {
        continue;
      }
      bytes = 0;

      for (const auto& it : futures) {
        it.wait();
      }
    }

    std::cerr << "[heron::] mapped sequences "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    j = i + 1;
  }

  timer.Stop();

  timer.Start();

  std::size_t zeros = 0, predictable = 0, total = 0, snps = 0, solid = 0;
  std::vector<std::vector<std::uint32_t>> rsnps(piles.size());
  std::size_t i = 0;
  for (const auto& it : piles) {
    total += it.size();
    std::size_t j = 0;
    for (const auto& jt : it) {
      std::vector<double> counts = {
          static_cast<double>(jt.a),
          static_cast<double>(jt.c),
          static_cast<double>(jt.g),
          static_cast<double>(jt.t)
      };
      double sum = std::accumulate(counts.begin(), counts.end(), 0);
      if (use_frequencies) {
        for (auto& kt : counts) {
          kt /= sum;
        }
      }
      if (sum == 0.) {
        ++zeros;
      } else if (sum > 7.) {
        ++predictable;
        std::size_t variants = 0;
        for (const auto& it : counts) {
          if (use_frequencies) {
            if (low < it && it < high) {
              ++variants;
            }
          } else {
            if (it > coverage) {
              ++variants;
            }
          }
        }
        if (variants > 1) {
          rsnps[i].emplace_back(j);
          ++snps;
        } else {
          ++solid;
        }
      }
      ++j;
    }
    ++i;
  }

  std::cerr << "[heron::] num bases = " << total << std::endl;
  std::cerr << "[heron::] num bad bases = " << zeros << std::endl;
  std::cerr << "[heron::] num predictable bases = " << predictable << std::endl;
  std::cerr << "[heron::] num solid = " << solid << std::endl;
  std::cerr << "[heron::] num snps = " << snps << std::endl;

  std::cerr << "[heron::] calculated statistics "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  for (std::uint32_t i = 0; i < rsnps.size(); ++i) {
    if (rsnps[i].empty()) {
      continue;
    }
    std::cout << i;
    for (const auto& jt : rsnps[i]) {
      std::cout << " " << jt;
    }
    std::cout << std::endl;
  }

  std::cerr << "[heron::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
