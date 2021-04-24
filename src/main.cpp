// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <cstdint>
#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "edlib.h"  // NOLINT
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
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
  std::uint32_t num_threads = 1;

  std::string optstr = "t:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 't': num_threads = atoi(optarg); break;
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

  auto sequences = sparser->Parse(-1);

  std::cerr << "[heron::] parsed " << sequences.size() << " sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  struct Pile {
    std::uint8_t a;
    std::uint8_t c;
    std::uint8_t g;
    std::uint8_t t;
  };

  std::vector<std::vector<Pile>> piles(sequences.size());
  for (const auto& it : sequences) {
    piles[it->id].resize(it->inflated_len);
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
  ram::MinimizerEngine minimizer_engine{thread_pool, 15U, 5U};

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

              EdlibAlignResult result = edlibAlign(
                  lhs.c_str(), lhs.size(),
                  rhs.c_str(), rhs.size(),
                  edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));  // NOLINT
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
                    case 1: ++lhs_pos; break;
                    case 2: ++rhs_pos; break;
                    default: break;
                  }
                }
              }
              edlibFreeAlignResult(result);
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
  for (const auto& it : piles) {
    total += it.size();
    for (const auto& jt : it) {
      std::vector<std::uint32_t> counts = { jt.a, jt.c, jt.g, jt.t };
      double sum = std::accumulate(counts.begin(), counts.end(), 0.);
      if (sum == 0.) {
        ++zeros;
      } else if (sum > 9.) {
        ++predictable;
        std::sort(counts.begin(), counts.end());
        if (counts[0] / sum > 0.5 && counts[1] / sum > 0.25) {
          ++snps;
        } else if (counts[0] / sum > 0.75) {
          ++solid;
        }
      }
    }
  }

  std::cerr << "[heron::] num bases = " << total << std::endl;
  std::cerr << "[heron::] num bad bases = " << zeros << std::endl;
  std::cerr << "[heron::] num predictable bases = " << predictable << std::endl;
  std::cerr << "[heron::] num solid = " << solid << std::endl;
  std::cerr << "[heron::] num snps = " << snps << std::endl;

  std::cerr << "[heron::] calculated statistics "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  std::cerr << "[heron::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
