//
//  csvwriter.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include "csvwriter.h"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <stdexcept>

#include <unistd.h>

std::string csvField(std::string_view s) {
  if (s.find_first_of(",\"\n") == std::string_view::npos)
    return std::string(s);
  std::string out = "\"";
  for (char c : s) {
    if (c == '"')
      out += '"';
    out += c;
  }
  out += '"';
  return out;
}

CsvWriter::CsvWriter(std::filesystem::path outputPath, std::string headerLine,
                     unsigned numThreads)
    : outputPath_(std::move(outputPath)), headerLine_(std::move(headerLine)),
      maxShards_(std::max<unsigned>(1, std::min(numThreads, MAX_SHARDS))) {
  static std::atomic<uint64_t> nextId{1};
  id_ = nextId.fetch_add(1);
}

void CsvWriter::writeRow(const std::string &row) {
  Shard &shard = shardForThisThread();
  std::lock_guard<std::mutex> lock(shard.mutex);
  shard.buffer += row;
  shard.buffer += '\n';
  if (shard.buffer.size() >= FLUSH_THRESHOLD)
    flush(shard);
}

void CsvWriter::finalize() {
  std::ofstream out(outputPath_, std::ios::trunc);
  if (!out)
    throw std::runtime_error("CsvWriter: failed to open " +
                             outputPath_.string() + " for writing");
  out << headerLine_ << "\n";
  for (auto &shard : shards_) {
    std::lock_guard<std::mutex> lock(shard->mutex);
    flush(*shard);
    shard->file.close();
    std::ifstream in(shard->path);
    if (!in)
      throw std::runtime_error("CsvWriter: failed to reopen " +
                               shard->path.string() + " for consolidation");
    out << in.rdbuf();
    in.close();
    std::filesystem::remove(shard->path);
  }
  if (!out)
    throw std::runtime_error("CsvWriter: failed while writing " +
                             outputPath_.string());
  std::cerr << "[+] Wrote CSV output to " << outputPath_.string() << "\n";
}

void CsvWriter::flush(Shard &shard) {
  if (shard.buffer.empty())
    return;
  shard.file << shard.buffer;
  if (!shard.file)
    throw std::runtime_error("CsvWriter: failed writing to " +
                             shard.path.string());
  shard.buffer.clear();
}

CsvWriter::Shard &CsvWriter::shardForThisThread() {
  static thread_local uint64_t cachedFor = 0;
  static thread_local Shard *cached = nullptr;
  if (cached && cachedFor == id_)
    return *cached;
  cachedFor = id_;
  std::lock_guard<std::mutex> lock(registryMutex_);
  if (shards_.size() < maxShards_) {
    auto shard = std::make_unique<Shard>();
    shard->path = outputPath_.string() + ".tmp." +
                 std::to_string(getpid()) + ".shard" +
                 std::to_string(shards_.size());
    shard->file.open(shard->path, std::ios::trunc);
    if (!shard->file)
      throw std::runtime_error("CsvWriter: failed to open " +
                               shard->path.string() + " for writing");
    shards_.push_back(std::move(shard));
  }
  cached = shards_[nextShard_++ % shards_.size()].get();
  return *cached;
}

// Minimal RFC-4180 field parser (quotes, doubled-quote escaping) -- needed
// to read our OWN --output file back on resume, since witness_pairsig/
// depends_on may have been written through csvField() and can contain
// commas.
std::vector<std::string> parseCsvLine(const std::string &line) {
  std::vector<std::string> fields;
  size_t i = 0;
  while (i <= line.size()) {
    std::string field;
    if (i < line.size() && line[i] == '"') {
      ++i;
      while (i < line.size()) {
        if (line[i] == '"') {
          if (i + 1 < line.size() && line[i + 1] == '"') {
            field += '"';
            i += 2;
          } else {
            ++i;
            break;
          }
        } else {
          field += line[i++];
        }
      }
    } else {
      while (i < line.size() && line[i] != ',')
        field += line[i++];
    }
    fields.push_back(field);
    if (i < line.size() && line[i] == ',') {
      ++i;
      continue;
    }
    break;
  }
  return fields;
}
