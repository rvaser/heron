// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "bioparser/fasta_parser.hpp"
#include "biosoup/sequence.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

static struct option options[] = {
    {"frequency", required_argument, nullptr, 'f'},
    {"seed", required_argument, nullptr, 's'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[heron::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz"
            << std::endl;
  return nullptr;
}


void Help() {
  std::cout <<
      "usage: mutate [options ...] <reference>\n"
      "\n"
      "  # default output is to stdout\n"
      "  <reference>\n"
      "    input file in FASTA format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -f, --frequency <double>\n"
      "      default: 0.01\n"
      "      mutation frequency\n"
      "    -s, --seed <int>\n"
      "      default: 42\n"
      "      random generator seed\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  double frequency = 0.01;
  int seed = 42;

  std::string optstr = "s:f:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'f': frequency = std::atof(optarg); break;
      case 's': seed = std::atoi(optarg); break;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc) {
    std::cerr << "[heron::Mutate] error: missing input file!" << std::endl;
    return 1;
  }

  auto rparser = CreateParser(argv[optind]);
  if (rparser == nullptr) {
    return 1;
  }
  std::vector<std::unique_ptr<biosoup::Sequence>> reference;
  try {
    reference = rparser->Parse(-1, false);
  } catch (const std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  std::mt19937 generator(seed);
  std::uniform_real_distribution<> distribution(0., 1.);

  for (const auto& it : reference) {
    std::cerr << it->id << " ";
    for (std::size_t j = 0; j < it->data.size(); ++j) {
      it->data[j] = std::toupper(it->data[j]);
      if (distribution(generator) <= frequency) {
        std::cerr << j << " ";
        switch (it->data[j]) {
          case 'A': it->data[j] = 'T'; break;
          case 'C': it->data[j] = 'G'; break;
          case 'G': it->data[j] = 'C'; break;
          case 'T': it->data[j] = 'A'; break;
          default: break;
        }
      }
    }
    std::cerr << std::endl;

    std::cout << ">" << it->name << " mutated" << std::endl
              << it->data << std::endl;
  }

  return 0;
}
