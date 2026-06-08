#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "base/core.h"
#include "base/rain.h"

using lms::core::ModelMeta;
using lms::rain::Station;
using lms::rain::Rainfall;
using lms::rain::RainfallEvent;

template <typename T, size_t kstation_cnt>
class RainBuilder {
 public:
  RainBuilder& from_directory(const std::string& dir, const std::string& file_name) {
    dir_ = dir;
    file_name_ = file_name;
    return *this;
  }

  void BuildAll(const ModelMeta<T>& meta) {
    BuildStationCsv();
    BuildRainCsv(meta.time_interval_s_);
    for (auto& station : stations_) {
      station.GeoPos2RasterPos(meta);
    }
  }

 private:
  std::filesystem::path dir_;
  std::string file_name_;
  std::vector<Station<T>> stations_;
  std::vector<RainfallEvent<T, kstation_cnt>> events_;

  void BuildStationCsv() {
    if (dir_.empty()) {
      throw std::runtime_error("StationBuilder: model directory not set");
    }
    std::ifstream file((dir_ / file_name_).string());
    if (!file.is_open()) {
      throw std::runtime_error("StationBuilder: Station data file not found");
    }

    std::string line;
    // read header (skip for now)
    std::getline(file, line);

    while (std::getline(file, line)) {
      std::stringstream ss(line);
      std::string id_str, name, geo_lat_str, geo_long_str;
      if (std::getline(ss, id_str, ',') && std::getline(ss, name, ',') &&
          std::getline(ss, geo_lat_str, ',') && std::getline(ss, geo_long_str, ',')) {
        std::size_t id = std::stoull(id_str);
        T geo_lat = static_cast<T>(std::stod(geo_lat_str));
        T geo_long = static_cast<T>(std::stod(geo_long_str));
        stations_.emplace_back(id, geo_lat, geo_long, name, ModelMeta<T>{});
      }
    }
  }

  void BuildRainCsv(std::size_t time_interval_s) {
    if (dir_.empty()) {
      throw std::runtime_error("RainBuilder: model directory not set");
    }
    std::ifstream file((dir_ / file_name_).string());
    if (!file.is_open()) {
      throw std::runtime_error("RainBuilder: rain data file not found");
    }

    std::string line;
    // read header: ID,Q,<station_name_1>,<station_name_2>,...
    std::getline(file, line);
    std::stringstream header_ss(line);
    std::string cell;
    std::vector<std::string> headers;
    while (std::getline(header_ss, cell, ',')) {
      headers.push_back(cell);
    }

    // Build mapping: CSV column index -> station index in stations_ array
    // headers[0] = "ID", headers[1] = "Q", headers[2..] = station names
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
      std::stringstream ss(line);
      std::string value;
      int col = 0;
      Rainfall<T, kstation_cnt> rainfall;
      rainfall.rainfall_row_.fill(T {0});

      while (std::getline(ss, value, ',')) {
        // col 0 = ID (time step), col 1 = Q (flow) -- skip
        // col 2.. = station rainfall values
        if (col >= 2 && col_to_station_idx[col] >= 0) {
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
