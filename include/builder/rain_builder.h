#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "base/core.h"
#include "base/rain.h"

namespace lms {
namespace rain {

using lms::core::ModelMeta;
using lms::rain::Station;
using lms::rain::Rainfall;
using lms::rain::RainfallEvent;

template <typename T, size_t kstation_cnt>
class RainBuilder {
 public:
  RainBuilder& from_directory(const std::string& dir) {
    dir_ = dir;
    return *this;
  }

  void BuildAll(const std::string& station_file, const std::string& rain_file, const ModelMeta<T>& meta) {
    BuildStationCsv(station_file);
    BuildRainCsv(rain_file, meta.time_interval_s_);
    for (auto& station : stations_) {
      station.GeoPos2RasterPos(meta);
    }
  }

  const std::vector<Station<T>>& stations() const { return stations_; }
  const std::vector<RainfallEvent<T, kstation_cnt>>& events() const { return events_; }

  std::vector<std::vector<T>> GetRainfallMatrix() const {
    std::vector<std::vector<T>> result;
    if (events_.empty()) return result;
    
    const auto& event = events_[0];
    result.reserve(event.rainfall_event_.size());
    for (const auto& rf : event.rainfall_event_) {
      std::vector<T> row(rf.rainfall_row_.begin(), rf.rainfall_row_.end());
      result.push_back(std::move(row));
    }
    return result;
  }

 private:
  std::filesystem::path dir_;
  std::vector<Station<T>> stations_;
  std::vector<RainfallEvent<T, kstation_cnt>> events_;

  void BuildStationCsv(const std::string& file_name) {
    std::ifstream file((dir_ / file_name).string());
    if (!file.is_open()) {
      throw std::runtime_error("RainBuilder: Station data file not found: " + file_name);
    }

    std::string line;
    // read header: ID,NAME,SID,X,Y,ENNAME
    std::getline(file, line);

    while (std::getline(file, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);
      std::vector<std::string> cols;
      std::string cell;
      while (std::getline(ss, cell, ',')) {
        cols.push_back(cell);
      }
      if (cols.size() >= 5) {
        std::size_t id = std::stoull(cols[0]);
        std::string name = cols[1]; // Use NAME for mapping
        T x = static_cast<T>(std::stod(cols[3]));
        T y = static_cast<T>(std::stod(cols[4]));
        stations_.emplace_back(id, y, x, name, ModelMeta<T>{}); // Assuming Y=lat, X=long
      }
    }
  }

  void BuildRainCsv(const std::string& file_name, std::size_t time_interval_s) {
    std::ifstream file((dir_ / file_name).string());
    if (!file.is_open()) {
      throw std::runtime_error("RainBuilder: rain data file not found: " + file_name);
    }

    std::string line;
    // read header: ID,Q,<station_name_1>,<station_name_2>,...
    if (!std::getline(file, line)) return;
    std::stringstream header_ss(line);
    std::string cell;
    std::vector<std::string> headers;
    while (std::getline(header_ss, cell, ',')) {
      headers.push_back(cell);
    }

    // Build mapping: CSV column index -> station index in stations_ array
    std::vector<int> col_to_station_idx(headers.size(), -1);
    for (std::size_t col = 2; col < headers.size(); ++col) {
      for (std::size_t si = 0; si < stations_.size(); ++si) {
        if (stations_[si].name_ == headers[col]) {
          col_to_station_idx[col] = static_cast<int>(si);
          break;
        }
      }
    }

    RainfallEvent<T, kstation_cnt> event;
    event.duration_ = 0;
    event.time_interval_s_ = time_interval_s;
    while (std::getline(file, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);
      std::string value;
      int col = 0;
      Rainfall<T, kstation_cnt> rainfall;
      rainfall.rainfall_row_.fill(T {0});

      while (std::getline(ss, value, ',')) {
        if (col >= 2 && col < (int)col_to_station_idx.size() && col_to_station_idx[col] >= 0) {
          rainfall.rainfall_row_[static_cast<std::size_t>(col_to_station_idx[col])] = static_cast<T>(std::stod(value));
        }
        ++col;
      }
      event.rainfall_event_.push_back(rainfall);
      ++event.duration_;
    }

    events_.push_back(std::move(event));
  }
};

} // namespace rain
} // namespace lms
