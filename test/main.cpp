// Driver to exercise the GDAL I/O layer against real data and print a summary,
// so raster input can be inspected without the (WIP) solver.
//
//   ./lms_io_test [model_dir]   (default: data/meizhou/model)

#include <algorithm>
#include <cstdio>
#include <limits>
#include <map>
#include <optional>
#include <string>

#include "base/model.h"
using T = double;
using lms::base::raster;
using lms::core::Label;
static const char* LabelName(Label l) {
  switch (l) {
    case Label::Soil:
      return "Soil";
    case Label::Channel:
      return "Channel";
    case Label::Reservoir:
      return "Reservoir";
  }
  return "?";
}

// Stats over active cells, optionally restricted to a single Label (used for
// channel-only parameters, which are NoData on soil cells).
template <typename Acc>
static void PrintFieldStats(const char* name, const ConstRaster<T>& cr, Acc acc,
                            std::optional<Label> only = std::nullopt) {
  double mn = std::numeric_limits<double>::max();
  double mx = std::numeric_limits<double>::lowest();
  double sum = 0.0;
  std::size_t cnt = 0;
  for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
    if (!cr.active[i]) continue;
    if (only && (*cr.cells)[i].label != *only) continue;
    const double val = static_cast<double>(acc((*cr.cells)[i]));
    mn = std::min(mn, val);
    mx = std::max(mx, val);
    sum += val;
    ++cnt;
  }
  if (cnt == 0) {
    std::printf("  %-10s : (no matching cells)\n", name);
    return;
  }
  std::printf("  %-10s : min=%-12.5g max=%-12.5g mean=%-12.5g\n", name, mn, mx, sum / cnt);
}

int main(int argc, char** argv) {
  const std::string dir = (argc > 1) ? argv[1] : "data/meizhou/model";

  try {
    // Raw label histogram: verify the label.tif -> Label mapping assumption.
    Raster<double> raw_label = ReadBand<double>(dir + "/label.tif");
    std::map<int, std::size_t> raw_hist;
    for (std::size_t i = 0; i < raw_label.data.size(); ++i) {
      if (raw_label.is_nodata_at(i)) continue;
      raw_hist[static_cast<int>(raw_label.data[i])]++;
    }

    ModelBuilder<T> builder;
    const ConstRaster<T> cr = builder.from_directory(dir).build_const_raster();

    std::size_t active = 0;
    std::map<int, std::size_t> label_counts;
    for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
      if (!cr.active[i]) continue;
      ++active;
      label_counts[static_cast<int>((*cr.cells)[i].label)]++;
    }

    std::printf("== Model I/O summary (%s) ==\n", dir.c_str());
    std::printf("grid        : %zu x %zu  (cell_size = %zu m)\n", cr.meta.width_, cr.meta.heigh_, cr.meta.cell_size_);
    std::printf("cells       : %zu total, %zu active (in basin)\n", (*cr.cells).size(), active);

    std::printf("raw label.tif values:\n");
    for (const auto& [val, c] : raw_hist) std::printf("  value %d : %zu cells\n", val, c);

    std::printf("decoded labels:\n");
    for (const auto& [val, c] : label_counts)
      std::printf("  %-10s : %zu cells\n", LabelName(static_cast<Label>(val)), c);

    std::printf("field stats over active cells:\n");
    PrintFieldStats("slop", cr, [](const ConstParam<T>& c) { return c.slop; });
    PrintFieldStats("ks", cr, [](const ConstParam<T>& c) { return c.ks; });
    PrintFieldStats("sat", cr, [](const ConstParam<T>& c) { return c.sat; });
    PrintFieldStats("fc", cr, [](const ConstParam<T>& c) { return c.fc; });
    PrintFieldStats("wl", cr, [](const ConstParam<T>& c) { return c.wl; });
    PrintFieldStats("zs", cr, [](const ConstParam<T>& c) { return c.zs; });
    PrintFieldStats("b", cr, [](const ConstParam<T>& c) { return c.b; });
    PrintFieldStats("n", cr, [](const ConstParam<T>& c) { return c.n; });
    PrintFieldStats("v", cr, [](const ConstParam<T>& c) { return c.v; });

    std::printf("channel-only params (over Channel cells):\n");
    PrintFieldStats("bs", cr, [](const ConstParam<T>& c) { return c.bs; }, Label::Channel);
    PrintFieldStats("bw", cr, [](const ConstParam<T>& c) { return c.bw; }, Label::Channel);
    PrintFieldStats("manning", cr, [](const ConstParam<T>& c) { return c.manning; }, Label::Channel);

    // Dump one active cell so the assembled ConstParam can be eyeballed.
    for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
      if (!cr.active[i]) continue;
      const ConstParam<T>& c = (*cr.cells)[i];
      std::printf(
          "sample active cell #%zu: label=%s d8=%d slop=%g ks=%g sat=%g fc=%g "
          "bw=%g manning=%g\n",
          i, LabelName(c.label), static_cast<int>(c.d8), static_cast<double>(c.slop), static_cast<double>(c.ks),
          static_cast<double>(c.sat), static_cast<double>(c.fc), static_cast<double>(c.bw),
          static_cast<double>(c.manning));
      break;
    }

    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
