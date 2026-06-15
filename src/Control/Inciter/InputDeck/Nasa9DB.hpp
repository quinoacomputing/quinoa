#pragma once
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

constexpr double NASA9_RU = 8.31446261815324; // J/(mol·K)

struct N9Interval {
  double Tlow  = 0.0;
  double Thigh = 0.0;
  std::array<double,9> a{};
};

struct N9Species {
  std::string name;
  double Mw         = 0.0;  // kg/mol
  double Hf298_mol  = 0.0;  // J/mol (formation)
  double Hf298_mass = 0.0;  // J/kg (formation)
  std::vector<N9Interval> intervals; // 3

  double R() const { return NASA9_RU / Mw; }   // J/(kg·K)

  // Backward-compatible accessors used elsewhere in Quinoa
  const std::vector<N9Interval>& getIntervals() const { return intervals; }

  std::size_t nIntervals() const { return intervals.size(); }

  const N9Interval& intervalByIndex(std::size_t i) const {
    return intervals.at(i);
  }

  const N9Interval& interval(double T) const {
    if (intervals.empty())
      Throw("No intervals for " + name);
    if (T <= intervals.front().Thigh) return intervals.front();
    if (T >= intervals.back().Tlow)   return intervals.back();
    for (const auto& I : intervals)
      if (T >= I.Tlow && T <= I.Thigh) return I;
    return intervals.back();
  }

  // Dimensionless NASA-9 forms
  static double Cp_over_R(double T, const N9Interval& I) {
    const double invT  = 1.0/T;
    const double invT2 = invT*invT;
    const double T2 = T*T, T3 = T2*T, T4 = T3*T;
    const auto& a = I.a;
    return a[0]*invT2 + a[1]*invT + a[2] + a[3]*T + a[4]*T2 + a[5]*T3 + a[6]*T4;
  }
  static double h_over_RT(double T, const N9Interval& I) {
    const double invT  = 1.0/T;
    const double invT2 = invT*invT;
    const double T2 = T*T, T3 = T2*T, T4 = T3*T;
    const auto& a = I.a;
    return -a[0]*invT2 + a[1]*std::log(T)*invT + a[2]
     + 0.5*a[3]*T + (a[4]*T2)/3.0 + (a[5]*T3)/4.0 + (a[6]*T4)/5.0
     + a[7]*invT;
  }
  static double s_over_R(double T, const N9Interval& I) {
    const double invT  = 1.0/T;
    const double invT2 = invT*invT;
    const double T2 = T*T, T3 = T2*T, T4 = T3*T;
    const auto& a = I.a;
    return -0.5*a[0]*invT2 - a[1]*invT + a[2]*std::log(T)
         + a[3]*T + 0.5*a[4]*T2 + (a[5]*T3)/3.0 + (a[6]*T4)/4.0
         + a[8];
  }
};

// Fixed-width 5×16 field reader (NASA9 style)
inline void n9_collect_fw(const std::string& line, std::vector<double>& out) {
  for (std::size_t f = 0; f < 5; ++f) {
    std::size_t start = f * 16;
    if (start >= line.size()) break;
    std::size_t end = std::min(start + 16, line.size());
    std::string chunk = line.substr(start, end - start);
    if (chunk.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    for (char& c : chunk) if (c=='D' || c=='d') c = 'E';
    char* e = nullptr;
    double v = std::strtod(chunk.c_str(), &e);
    if (e != chunk.c_str())
      out.push_back(v);
  }
}

// Whitespace float reader
inline void n9_collect_ws(const std::string& s, std::vector<double>& out) {
  std::istringstream ss(s);
  std::string tok;
  while (ss >> tok) {
    for (char& c : tok) if (c=='D' || c=='d') c='E';
    char* e = nullptr;
    double v = std::strtod(tok.c_str(), &e);
    if (e != tok.c_str())
      out.push_back(v);
  }
}

// Normalize to exactly 3 intervals
inline void n9_normalize_to_three_intervals(N9Species& sp) {
  if (sp.intervals.size() == 3) return;
  if (sp.intervals.empty()) Throw("No intervals to normalize for " + sp.name);

  if (sp.intervals.size() == 2) {
    const N9Interval last = sp.intervals.back();
    const double b = last.Tlow;
    const double c = last.Thigh;
    if (!(c > b)) Throw("Invalid bounds for " + sp.name + ": Tlow >= Thigh");
    const double mid = 0.5*(b + c);

    N9Interval I1 = last; I1.Tlow = b;   I1.Thigh = mid;
    N9Interval I2 = last; I2.Tlow = mid; I2.Thigh = c;

    sp.intervals.pop_back();
    sp.intervals.push_back(I1);
    sp.intervals.push_back(I2);
    return;
  }

  if (sp.intervals.size() == 1) {
    const N9Interval only = sp.intervals.front();
    const double a = only.Tlow;
    const double b = only.Thigh;
    if (!(b > a)) Throw("Invalid bounds for " + sp.name + ": Tlow >= Thigh");

    const double d  = (b - a) / 3.0;
    const double t1 = a + d;
    const double t2 = a + 2.0*d;

    N9Interval I0 = only; I0.Tlow = a;  I0.Thigh = t1;
    N9Interval I1 = only; I1.Tlow = t1; I1.Thigh = t2;
    N9Interval I2 = only; I2.Tlow = t2; I2.Thigh = b;

    sp.intervals.clear();
    sp.intervals.push_back(I0);
    sp.intervals.push_back(I1);
    sp.intervals.push_back(I2);
    return;
  }

  Throw("Unsupported number of intervals for " + sp.name + ": " +
        std::to_string(sp.intervals.size()));
}

// NASA9 reader (now supports 1–3 intervals, normalized to 3)
inline N9Species read_nasa9_species(const std::string& file,
                                   const std::string& targetName)
{
  std::ifstream in(file);
  if (!in) Throw("Cannot open nasa9 file: " + file);

  std::string line;
  while (std::getline(in, line)) {
    std::string nameLine = line;
    auto p1 = nameLine.find_first_not_of(" \t\r\n");
    if (p1 == std::string::npos) continue;

    std::istringstream ns(nameLine.substr(p1));
    std::string spName;
    ns >> spName;
    if (spName != targetName) continue;

    std::string comp;
    if (!std::getline(in, comp))
      Throw("EOF after header for " + targetName);

    std::vector<double> nums;
    n9_collect_ws(comp, nums);
    if (nums.size() < 2)
      Throw("Cannot parse Mw/Hf line for " + targetName);

    double Hf298   = nums.back();            // J/mol
    double Mw_gmol = nums[nums.size()-2];    // g/mol
    double Mw      = Mw_gmol * 1e-3;         // kg/mol

    N9Species sp;
    sp.name        = targetName;
    sp.Mw          = Mw;
    sp.Hf298_mol   = Hf298;
    sp.Hf298_mass  = Hf298 / Mw;

    const std::size_t nIntervals = static_cast<std::size_t>(std::lround(nums[0]));
    if (nIntervals < 1 || nIntervals > 3)
      Throw("Only species with 1, 2, or 3 temperature intervals are supported. " +
            targetName + " has " + std::to_string(nIntervals));

    for (std::size_t iv = 0; iv < nIntervals; ++iv) {
      std::string hdr;
      while (true) {
        if (!std::getline(in, hdr))
          Throw("EOF in header for " + targetName);
        std::string tmp = hdr;
        auto q1 = tmp.find_first_not_of(" \t\r\n");
        if (q1 == std::string::npos || tmp[q1]=='!') continue;
        break;
      }

      std::vector<double> hnums;
      n9_collect_ws(hdr, hnums);
      if (hnums.size() < 2)
        Throw("Bad interval header for " + targetName);
      double Tmin = hnums[0];
      double Tmax = hnums[1];

      std::string c1, c2;
      if (!std::getline(in, c1) || !std::getline(in, c2))
        Throw("EOF reading coeffs for " + targetName);

      std::vector<double> coeffs;
      coeffs.reserve(10);
      n9_collect_fw(c1, coeffs);
      n9_collect_fw(c2, coeffs);
      if (coeffs.size() != 9 && coeffs.size() != 10)
        Throw("Expected 9 or 10 coeffs for " + targetName +
              ", got " + std::to_string(coeffs.size()));

      N9Interval I;
      I.Tlow  = Tmin;
      I.Thigh = Tmax;

      // Mutation++/NASA9 files may appear in two layouts:
      //   9 coeffs:  a0 ... a6, a7(H), a8(S)
      //  10 coeffs: a0 ... a6, dummy, a7(H), a8(S)
      // The dummy field is not a NASA9 thermodynamic coefficient.
      for (std::size_t k = 0; k < 7; ++k) I.a[k] = coeffs[k];
      if (coeffs.size() == 10) {
        I.a[7] = coeffs[8];
        I.a[8] = coeffs[9];
      } else {
        I.a[7] = coeffs[7];
        I.a[8] = coeffs[8];
      }

      sp.intervals.push_back(I);
    }

    n9_normalize_to_three_intervals(sp);
    return sp;
  }

  Throw("Species " + targetName + " not found in " + file);
}

inline bool try_read_nasa9_species(const std::string& file,
                                  const std::string& targetName,
                                  N9Species& out)
{
  std::ifstream in(file);
  if (!in) Throw("Cannot open nasa9 file: " + file);

  std::string line;
  while (std::getline(in, line)) {
    std::string nameLine = line;
    auto p1 = nameLine.find_first_not_of(" \t\r\n");
    if (p1 == std::string::npos) continue;

    std::istringstream ns(nameLine.substr(p1));
    std::string spName;
    ns >> spName;
    if (spName != targetName) continue;

    // ---- Found species ----
    std::string comp;
    if (!std::getline(in, comp))
      Throw("EOF after header for " + targetName);

    std::vector<double> nums;
    n9_collect_ws(comp, nums);
    if (nums.size() < 2)
      Throw("Cannot parse Mw/Hf line for " + targetName);

    double Hf298   = nums.back();            // J/mol
    double Mw_gmol = nums[nums.size()-2];    // g/mol
    double Mw      = Mw_gmol * 1e-3;         // kg/mol

    N9Species sp;
    sp.name        = targetName;
    sp.Mw          = Mw;
    sp.Hf298_mol   = Hf298;
    sp.Hf298_mass  = Hf298 / Mw;

    const std::size_t nIntervals = static_cast<std::size_t>(std::lround(nums[0]));
    if (nIntervals < 1 || nIntervals > 3)
      Throw("Only species with 1, 2, or 3 temperature intervals are supported. " +
            targetName + " has " + std::to_string(nIntervals));

    for (std::size_t iv = 0; iv < nIntervals; ++iv) {
      std::string hdr;
      while (true) {
        if (!std::getline(in, hdr))
          Throw("EOF in header for " + targetName);
        std::string tmp = hdr;
        auto q1 = tmp.find_first_not_of(" \t\r\n");
        if (q1 == std::string::npos || tmp[q1]=='!') continue;
        break;
      }

      std::vector<double> hnums;
      n9_collect_ws(hdr, hnums);
      if (hnums.size() < 2)
        Throw("Bad interval header for " + targetName);
      double Tmin = hnums[0];
      double Tmax = hnums[1];

      std::string c1, c2;
      if (!std::getline(in, c1) || !std::getline(in, c2))
        Throw("EOF reading coeffs for " + targetName);

      std::vector<double> coeffs;
      coeffs.reserve(10);
      n9_collect_fw(c1, coeffs);
      n9_collect_fw(c2, coeffs);
      if (coeffs.size() != 9 && coeffs.size() != 10)
        Throw("Expected 9 or 10 coeffs for " + targetName +
              ", got " + std::to_string(coeffs.size()));

      N9Interval I;
      I.Tlow  = Tmin;
      I.Thigh = Tmax;

      // Mutation++/NASA9 files may appear in two layouts:
      //   9 coeffs:  a0 ... a6, a7(H), a8(S)
      //  10 coeffs: a0 ... a6, dummy, a7(H), a8(S)
      // The dummy field is not a NASA9 thermodynamic coefficient.
      for (std::size_t k = 0; k < 7; ++k) I.a[k] = coeffs[k];
      if (coeffs.size() == 10) {
        I.a[7] = coeffs[8];
        I.a[8] = coeffs[9];
      } else {
        I.a[7] = coeffs[7];
        I.a[8] = coeffs[8];
      }

      sp.intervals.push_back(I);
    }

    n9_normalize_to_three_intervals(sp);
    out = std::move(sp);
    return true;
  }

  return false;
}
