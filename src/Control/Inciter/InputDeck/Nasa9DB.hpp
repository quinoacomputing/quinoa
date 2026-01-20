#pragma once
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cmath>
#include <stdexcept>

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
  std::vector<N9Interval> intervals; // should be 3

  double R() const { return NASA9_RU / Mw; }   // J/(kg·K)

  const std::vector<N9Interval>& getIntervals() const { return intervals; }

  std::size_t nIntervals() const { return intervals.size(); }

  const N9Interval& intervalByIndex(std::size_t i) const {
    return intervals.at(i);
  }

  const N9Interval& interval(double T) const {
    if (intervals.empty())
      throw std::runtime_error("No intervals for " + name);
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
    return -a[0]*invT2 + a[1]*std::log(T) + a[2]
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

  // Mass-based properties
  double Cp(double T) const { // J/(kg·K)
    const auto& I = interval(T);
    return Cp_over_R(T,I) * NASA9_RU / Mw;
  }
  double h(double T) const {  // J/kg (sensible)
    const auto& I = interval(T);
    return h_over_RT(T,I) * NASA9_RU * T / Mw;
  }
  double s(double T) const {  // J/(kg·K)
    const auto& I = interval(T);
    return s_over_R(T,I) * NASA9_RU / Mw;
  }
  double h_ref_298() const { return h(298.15); } // J/kg
};

// Fixed-width 5×16 field reader (NASA style)
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

// Whitespace float reader (for comp/header lines)
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

// Read ONLY ONE SPECIES block from nasa9.dat
inline N9Species read_nasa9_species(const std::string& file,
                                    const std::string& targetName)
{
  std::ifstream in(file);
  if (!in)
    throw std::runtime_error("Cannot open nasa9 file: " + file);

  std::string line;
  while (std::getline(in, line)) {
    // strip leading spaces
    std::string nameLine = line;
    auto p1 = nameLine.find_first_not_of(" \t\r\n");
    if (p1 == std::string::npos) continue;
    // species name is first token
    std::istringstream ns(nameLine.substr(p1));
    std::string spName;
    ns >> spName;
    if (spName != targetName) continue; // not the one we want

    // ---- Found the species header ----

    // comp/meta line
    std::string comp;
    if (!std::getline(in, comp))
      throw std::runtime_error("EOF after header for " + targetName);

    std::vector<double> nums;
    n9_collect_ws(comp, nums);
    if (nums.size() < 2)
      throw std::runtime_error("Cannot parse Mw/Hf line for " + targetName);

    double Hf298   = nums.back();            // J/mol
    double Mw_gmol = nums[nums.size()-2];    // g/mol
    double Mw      = Mw_gmol * 1e-3;         // kg/mol

    N9Species sp;
    sp.name        = targetName;
    sp.Mw          = Mw;
    sp.Hf298_mol   = Hf298;
    sp.Hf298_mass  = Hf298 / Mw;

    // 3 intervals: header + 2 coeff lines each
    for (std::size_t iv = 0; iv < 3; ++iv) {
      std::string hdr;
      // skip comments/blank lines between blocks
      while (true) {
        if (!std::getline(in, hdr))
          throw std::runtime_error("EOF in header for " + targetName);
        std::string tmp = hdr;
        auto q1 = tmp.find_first_not_of(" \t\r\n");
        if (q1 == std::string::npos || tmp[q1]=='!') continue;
        break;
      }

      std::string thdr = hdr;
      // parse Tmin Tmax from header
      std::vector<double> hnums;
      n9_collect_ws(thdr, hnums);
      if (hnums.size() < 2)
        throw std::runtime_error("Bad interval header for " + targetName);
      double Tmin = hnums[0];
      double Tmax = hnums[1];

      // coeff lines
      std::string c1, c2;
      if (!std::getline(in, c1) || !std::getline(in, c2))
        throw std::runtime_error("EOF reading coeffs for " + targetName);

      std::vector<double> coeffs;
      coeffs.reserve(10);
      n9_collect_fw(c1, coeffs);
      n9_collect_fw(c2, coeffs);
      if (coeffs.size() < 9)
        throw std::runtime_error("Fewer than 9 coeffs for " + targetName);

      N9Interval I;
      I.Tlow  = Tmin;
      I.Thigh = Tmax;
      for (std::size_t k = 0; k < 9; ++k)
        I.a[k] = coeffs[k];

      sp.intervals.push_back(I);
    }

    return sp;
  }

  throw std::runtime_error("Species " + targetName + " not found in " + file);
}
