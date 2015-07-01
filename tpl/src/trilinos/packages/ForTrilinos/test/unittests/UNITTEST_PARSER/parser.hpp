#ifndef PARSER_HH
#define PARSER_HH

#include <vector>
#include <string>
#include <map>

class Parser
{
  public:
    Parser(const std::string &prefix)
      : _prefix(prefix) { }
    virtual ~Parser() { }

    virtual void init() = 0;
    std::vector<std::string> process_line(const std::string &line);

  protected:
    std::string _prefix;
    typedef void (*proc_fun_t)(const std::string &line, const std::string &term,
        const std::vector<std::string> &args, std::vector<std::string> &result);
    void add_match(const std::string &term, proc_fun_t fun);
    std::string strip_space(const std::string &par);
    bool is_preproc(const std::string &line);
    bool is_blank(const std::string &line);

  private:
    typedef std::map<std::string, proc_fun_t> proc_fun_map_t;

    static proc_fun_map_t &getMap();
    bool is_match(const std::string &line, std::string &term,
        std::vector<std::string> &args);
    bool extract_match(const std::string &line, std::string &match,
        std::string &args);
    void break_args(const std::string &argstr, std::vector<std::string> &args);
    void proc_match(const std::string &line, const std::string &term,
        std::vector<std::string> &args, std::vector<std::string> &result);
    virtual void proc_plain(const std::string &line, std::vector<std::string> &result);
};

class ImplParser : public Parser
{
  public:
    ImplParser(const std::string &prefix)
      : Parser(prefix) { }
    virtual ~ImplParser() { }

    void init();
};

class CallsParser : public Parser
{
  public:
    CallsParser(const std::string &prefix)
      : Parser(prefix) { }
    virtual ~CallsParser() { }

    void init();

  private:
    void proc_plain(const std::string &line, std::vector<std::string> &result);
};

class ListParser : public Parser
{
  public:
    ListParser(const std::string &prefix)
      : Parser(prefix) { }
    virtual ~ListParser() { }

    void init();

  private:
    void proc_plain(const std::string &line, std::vector<std::string> &result);
};

class MpiListParser : public Parser
{
  public:
    MpiListParser(const std::string &prefix)
      : Parser(prefix) { }
    virtual ~MpiListParser() { }

    void init();

  private:
    void proc_plain(const std::string &line, std::vector<std::string> &result);
};

#endif /* PARSER_HH */
