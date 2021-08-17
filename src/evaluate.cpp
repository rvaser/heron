// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

static struct option options[] = {
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

void Help() {
  std::cout <<
      "usage: evaluate [options ...] <annotations> <annotations>\n"
      "\n"
      "  # default output is to stdout\n"
      "  <annotations>\n"
      "    file containing sequence variant calls separated with spaces,\n"
      "    starting with the sequence id\n"
      "\n"
      "  options:\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input_paths;

  std::string optstr = "h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
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

  if (input_paths.size() < 2) {
    std::cerr << "[heron::Evaluate] error: missing input file(s)!" << std::endl;
    return 1;
  }

  auto a = ParseAnnotations(input_paths[0]);
  if (a.empty()) {
    return 1;
  }

  auto b = ParseAnnotations(input_paths[1]);
  if (b.empty()) {
    return 1;
  }

  std::uint32_t max_id = 0;
  for (const auto& it : a) {
    max_id = std::max(max_id, it.first);
  }
  for (const auto& it : b) {
    max_id = std::max(max_id, it.first);
  }

  std::vector<double> jaccard;
  for (std::uint32_t i = 0; i < max_id + 1; ++i) {
    if (a.find(i) == a.end() && b.find(i) == b.end()) {
      continue;
    }
    const auto& it = a[i];
    const auto& jt = b[i];
    std::size_t intersection = 0;
    for (auto kt : it) {
      if (jt.find(kt) != jt.end()) {
        ++intersection;
      }
    }
    jaccard.emplace_back(
        intersection /
        static_cast<double>(it.size() + jt.size() - intersection));
  }
  std::sort(jaccard.begin(), jaccard.end());

  std::cout << jaccard.front() << " "
            << jaccard[jaccard.size() / 2] << " "
            << std::accumulate(jaccard.begin(), jaccard.end(), 0.) / jaccard.size() << " " // NOLINT
            << jaccard.back() << std::endl;

  return 0;
}
