// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "edlib.h"  // NOLINT
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unordered_map<std::uint32_t, std::unordered_set<std::uint32_t>> ParseAnnotations(  // NOLINT
    const std::string& path) {
  std::unordered_map<std::uint32_t, std::unordered_set<std::uint32_t>> dst;

  std::ifstream ifs(path);
  if (ifs.fail()) {
    std::cerr << "[heron::ParseAnnotations] error: unable to open file " << path
              << std::endl;
    return dst;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty()) {
      break;
    }
    std::istringstream iss(line);
    std::size_t id, pos;
    iss >> id;
    while (iss >> pos) {
      dst[id].emplace(pos);
    }
  }
  ifs.close();

  return dst;
}

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
      "usage: annotate [options ...] <sequences> <reference> <annotations>\n"
      "\n"
      "  # default output is to stdout\n"
      "  <sequences>/<reference>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "  <annotations>\n"
      "    file containing reference variant calls separated with spaces,\n"
      "    starting with the reference id\n"
      "\n"
      "  options:\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::size_t num_threads = 1;

  std::vector<std::string> input_paths;

  std::string optstr = "t:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 't': num_threads = std::atoi(optarg); break;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  for (std::int32_t i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.size() < 3) {
    std::cerr << "[heron::Annotate] error: missing input file(s)!" << std::endl;
    return 1;
  }

  auto sparser = CreateParser(input_paths[0]);
  if (sparser == nullptr) {
    return 1;
  }
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  try {
    sequences = sparser->Parse(-1);
  } catch (const std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  biosoup::NucleicAcid::num_objects = 0;

  auto rparser = CreateParser(input_paths[1]);
  if (rparser == nullptr) {
    return 1;
  }
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> reference;
  try {
    reference = rparser->Parse(-1);
  } catch (const std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  auto annotations = ParseAnnotations(input_paths[2]);
  if (annotations.empty()) {
    return 1;
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  ram::MinimizerEngine minimizer_engine{thread_pool};
  minimizer_engine.Minimize(reference.begin(), reference.end());
  minimizer_engine.Filter(0.001);

  auto overlap_length = [] (const biosoup::Overlap& o) -> std::uint32_t {
    return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
  };
  auto overlap_score = [&] (const biosoup::Overlap& o) -> double {
    return ((2. * overlap_length(o)) * o.score) / (overlap_length(o) + o.score);
  };

  std::vector<std::future<std::vector<std::uint32_t>>> futures;
  for (const auto& it : sequences) {
    futures.emplace_back(thread_pool->Submit(
        [&] (decltype(it) sequence) -> std::vector<std::uint32_t> {
          std::vector<std::uint32_t> dst;

          auto ovl = minimizer_engine.Map(sequence, false, false);
          if (ovl.empty()) {
            return dst;
          }

          std::sort(ovl.begin(), ovl.end(),
              [&] (const biosoup::Overlap& lhs,
                   const biosoup::Overlap& rhs) -> bool {
                return overlap_score(lhs) > overlap_score(rhs);
              });

          auto it = ovl.front();

          biosoup::NucleicAcid lhs_{"", sequence->InflateData(
              it.lhs_begin,
              it.lhs_end - it.lhs_begin)};
          if (!it.strand) {
            lhs_.ReverseAndComplement();
            auto tmp = it.lhs_begin;
            it.lhs_begin = sequence->inflated_len - it.lhs_end;
            it.lhs_end = sequence->inflated_len - tmp;
          }
          auto lhs = lhs_.InflateData();

          auto rhs = reference[it.rhs_id]->InflateData(
              it.rhs_begin,
              it.rhs_end - it.rhs_begin);

          auto result = edlibAlign(
              lhs.c_str(), lhs.size(),
              rhs.c_str(), rhs.size(),
              edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));  // NOLINT

          if (result.status == EDLIB_STATUS_OK) {
            dst.emplace_back(sequence->id);

            std::uint32_t lhs_pos = it.lhs_begin;
            std::uint32_t rhs_pos = it.rhs_begin;
            for (int i = 0; i < result.alignmentLength; ++i) {
              switch (result.alignment[i]) {
                case 0:
                case 3:
                  if (annotations[it.rhs_id].find(rhs_pos) !=
                      annotations[it.rhs_id].end()) {
                    dst.emplace_back(it.strand ?
                        lhs_pos :
                        sequence->inflated_len - lhs_pos - 1);
                  }
                  ++lhs_pos;
                  ++rhs_pos;
                  break;
                case 1: ++lhs_pos; break;
                case 2: ++rhs_pos; break;
                default: break;
              }
            }
          }
          edlibFreeAlignResult(result);

          return dst;
        },
        std::ref(it)));
  }
  for (auto& it : futures) {
    auto variants = it.get();
    if (variants.empty()) {
      continue;
    }

    std::sort(variants.begin() + 1, variants.end());
    for (const auto& jt : variants) {
      std::cout << jt << " ";
    } std::cout << std::endl;
  }

  return 0;
}
