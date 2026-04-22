// include/Plotting.h
#ifndef PLOTTING_H
#define PLOTTING_H

#include <string>
#include <vector>
#include <tuple>

/**
 * Draw a comparison plot with multiple model curves (lines) and optional experimental data points.
 * The plot is saved as both PDF and PNG with the given filename (without extension).
 *
 * @param filename  Base name of the output file (e.g., "proton_mt_distribution")
 * @param title     Title of the plot (appears at the top)
 * @param x_title   Label for the x‑axis
 * @param y_title   Label for the y‑axis
 * @param series    Vector of tuples defining model curves.
 *                  Each tuple: (x_data, y_data, line_color, line_style, legend_label)
 *                  The vectors must outlive the drawing call (they are only referenced).
 * @param exp_data  Vector of tuples defining experimental data points.
 *                  Each tuple: (x_data, y_data, marker_color, marker_style, legend_label)
 *                  The vectors must outlive the drawing call.
 * @param logy      If true, set y‑axis to logarithmic scale.
 * @param x_min     Minimum x‑axis limit (if > x_max, automatic scaling is used)
 * @param x_max     Maximum x‑axis limit
 * @param y_min     Minimum y‑axis limit (if > y_max, automatic scaling is used)
 * @param y_max     Maximum y‑axis limit
 */
void DrawComparisonPlot(
    const std::string& filename,
    const std::string& title,
    const std::string& x_title,
    const std::string& y_title,
    const std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>>& series,
    const std::vector<std::tuple<std::vector<double>*, std::vector<double>*, int, int, std::string>>& exp_data,
    bool logy = false,
    double x_min = 0, double x_max = 0,
    double y_min = 0, double y_max = 0
);

#endif // PLOTTING_H