#include <iostream>
#include "parser.hpp"
#include "exception.hpp"

std::vector<std::string> Parser::process_line(const std::string &line)
{
  std::string term;
  std::vector<std::string> args;
  std::vector<std::string> result;
  bool match = is_match(line, term, args);
  if (match) {
    proc_match(line, term, args, result);
  } else {
    proc_plain(line, result);
  }
  return result;
}

Parser::proc_fun_map_t &Parser::getMap()
{
  static proc_fun_map_t map;
  return map;
}

void Parser::add_match(const std::string &term, proc_fun_t fun)
{
  proc_fun_map_t &m = getMap();
  m[term] = fun;
}

bool Parser::is_match(const std::string &line, std::string &term,
    std::vector<std::string> &args)
{
  std::string my_term;
  std::string argstr;
  bool match = extract_match(line, my_term, argstr);
  if (match) {
    proc_fun_map_t &m = getMap();
    if (m.find(my_term) == m.end()) {
      return false;
    } else {
      term = my_term;
      break_args(argstr, args);
    }
  }
  return match;
}

bool Parser::extract_match(const std::string &line, std::string &match,
    std::string &args)
{
  size_t loc_prefix = line.find(_prefix);
  if (loc_prefix == line.npos)
    return false;

  size_t tmp = line.find("!");
  if ((tmp != line.npos) && (tmp < loc_prefix))
    return false;

  size_t loc_term = loc_prefix + _prefix.size()+1;
  size_t loc_lpar = line.find("(", loc_term);
  size_t loc_rpar = line.rfind(")");
  size_t loc_space = line.find(" ", loc_term);

  size_t loc_after_term;
  if ((loc_lpar != line.npos) && (loc_space != line.npos))
    loc_after_term = (loc_lpar < loc_space ? loc_lpar : loc_space);
  else if (loc_lpar != line.npos)
    loc_after_term = loc_lpar;
  else // even if it's npos
    loc_after_term = loc_space;
    
  std::string my_term(line, loc_term, loc_after_term-loc_term);
  match = my_term;

  if (loc_lpar != line.npos) {
    if (loc_rpar == line.npos) {
      throw ParserException("Syntax error in argument list");
    } else {
      std::string my_args(line, loc_lpar+1, loc_rpar-loc_lpar-1);
      args = my_args;
    }
  } else {
    std::string my_args("");
    args = my_args;
  }

  return true;
}

std::string Parser::strip_space(const std::string &par)
{
  size_t len = par.size();
  size_t firstgood = 0;
  size_t lastgood = len-1;
  if (len > 0) {
    for (size_t i=0; i<len; i++) {
      if (par[i] != ' ') break;
      firstgood = i+1;
    }
    for (size_t i=len-1; i>=0; i--) {
      if (par[i] != ' ') break;
      lastgood = i-1;
    }
    if ((firstgood != 0) || (lastgood != (len-1))) {
      std::string tmp(par, firstgood, lastgood-firstgood+1);
      return tmp;
    }
  }
  return par;
}

void Parser::break_args(const std::string &argstr, std::vector<std::string> &args)
{
  args.clear();
  if (argstr.size() == 0) return;
  std::string remaining = argstr;
  bool done = false;
  while (!done) {
    size_t loc_comma = remaining.find(",");
    size_t loc_lpar = remaining.find("(");
    if (loc_comma == remaining.npos) {
      args.push_back(strip_space(remaining));
      done = true;
    } else {
      if (loc_lpar != remaining.npos) {
        size_t loc_rpar = remaining.find(")", loc_lpar);
        loc_comma = remaining.find(",", loc_rpar);
        if (loc_comma == remaining.npos) {
          args.push_back(strip_space(remaining));
          done = true;
        } else {
          std::string first(remaining, 0, loc_comma);
          args.push_back(strip_space(first));
          std::string leftover(remaining, loc_comma+1);
          remaining = leftover;
        }
      } else {
        std::string first(remaining, 0, loc_comma);
        args.push_back(strip_space(first));
        std::string leftover(remaining, loc_comma+1);
        remaining = leftover;
      }
    }
  }
}

void Parser::proc_match(const std::string &line, const std::string &term,
    std::vector<std::string> &args, std::vector<std::string> &result)
{
  proc_fun_map_t &m = getMap();
  if (m.find(term) == m.end()) {
    throw ParserException("Internal error");
  } else {
    (m[term])(line, term, args, result);
  }
}

bool Parser::is_preproc(const std::string &line)
{
  std::string mline = strip_space(line);
  if (mline.size() > 0) {
    return (mline[0] == '#');
  } else {
    return false;
  }
}

bool Parser::is_blank(const std::string &line)
{
  std::string mline = strip_space(line);
  return (mline.size() == 0);
}

void Parser::proc_plain(const std::string &line, std::vector<std::string> &result)
{
  if (is_preproc(line)) {
    result.push_back(strip_space(line));
  } else {
    result.push_back(line);
  }
}

