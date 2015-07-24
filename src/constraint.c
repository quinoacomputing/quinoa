/*
 o---------------------------------------------------------------------o
 |
 | Ndiff
 |
 | Copyright (c) 2012+ laurent.deniau@cern.ch
 | Gnu General Public License
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     create constraints content
     print, scan constraints from file
 
 o---------------------------------------------------------------------o
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "args.h"
#include "utils.h"
#include "error.h"
#include "register.h"
#include "constraint.h"

#define T struct constraint
#define S struct slice

// ----- private

static void
printSlc(const S *s, FILE *out)
{
  if (slice_isFull(s)) {
    putc('*', out);
    return;
  }

  fprintf(out, "%u", slice_first(s));

  if (slice_isUnit(s)) return;

  if (slice_isInfinite(s)) {
    putc('-', out);
    putc('$', out);
  } else
   fprintf(out, "-%u", slice_last(s));

  if (slice_stride(s) != 1)
    fprintf(out, "/%u", slice_stride(s));    
}

static int
readSlcOrRng(S *s, FILE *in)
{
  int c, r = 1;
  uint first=0, last=0, stride=1;

  // skip spaces
  while((c = getc(in)) != EOF && isblank(c)) ;
  if (c == EOF) return EOF;

  // ('*'|num)
  if (c == '*') { last = UINT_MAX; goto finish; }
  else {
    ungetc(c, in);
    if (fscanf(in, "%u", &first) != 1) return EOF;
  }

  // (':'|'-')?
  c = getc(in);
       if (c == ':') r = 0;  // slice
  else if (c == '-') ;       // range
  else { ungetc(c, in); last = first; goto finish; }

  // ('$'|num)
  c = getc(in);
  if (c == '$') last = UINT_MAX;
  else {
    ungetc(c, in);
    if (fscanf(in, "%u", &last) != 1) return EOF;
  }

  // ('/'num)? 
  c = getc(in);
  if (c != '/') { ungetc(c, in); stride = 1; goto finish; }
  else
    if (fscanf(in, "%u", &stride) != 1) return EOF;

finish:
  if (r)
    *s = slice_initLastStride(first, last, stride);
  else
    *s = slice_initSizeStride(first, last, stride);

  trace("<-readSlcOrRng %u%c%u/%u", first, r ? '-' : ':', last, stride);

  return 0;
}

#define EPS_INVALID \
{ \
  cmd = eps_invalid; \
  trace("[%d] invalid command str='%s', buf='%s'", row, str, buf); \
  break; \
}

static int
readEps(struct eps *e, FILE *in, int row)
{
  int c = 0, n = 0, cmd = eps_invalid;
  char str[16], buf[MkConcat(1,MAXTAGLEN)], *end;

  while (1) {
    // parse next command
    *str = *buf = 0;
    n = fscanf(in, "%*[ \t]%16[^= \t\n\r!#]", str);
    str[sizeof str-1]=0;

    if (n == EOF || *str == 0) break;

// commands without '='
    if ((c = getc(in)) != '=') {
      ungetc(c, in);

           if (strcmp(str, "skip") == 0) {
        cmd |= eps_skip;  trace("[%d] skip", row);
      }
      else if (strcmp(str, "ign") == 0) {
        cmd |= eps_ign;   trace("[%d] ign", row);
      }
      else if (strcmp(str, "istr") == 0) {
        cmd |= eps_istr;  trace("[%d] istr", row);
      }
      else if (strcmp(str, "equ") == 0) {
        cmd |= eps_equ;   trace("[%d] equ", row);
      }
      else if (strcmp(str, "any") == 0) {
        cmd |= eps_any;   trace("[%d] any", row);
      }
      else if (strcmp(str, "all") == 0) {
        cmd &= ~eps_any;  trace("[%d] all", row);
      }
      else if (strcmp(str, "alt") == 0) {
        cmd |= eps_alt;  trace("[%d] alt", row);
      }
      else if (strcmp(str, "eval") == 0) {
        cmd |= eps_eval; trace("[%d] eval", row);
      }
      else if (strcmp(str, "nofail") == 0) {
        cmd |= eps_nofail; trace("[%d] nofail", row);
      }
      else if (strcmp(str, "large") == 0) {
        cmd |= eps_large; trace("[%d] large", row);
      }
      else if (strcmp(str, "small") == 0) {
        cmd &= ~eps_large; trace("[%d] small", row);
      }
      else if (strcmp(str, "trace") == 0) {
        cmd |= eps_trace; trace("[%d] trace", row);
      }
      else if (strcmp(str, "traceR") == 0) {
        cmd |= eps_trace | eps_traceR; trace("[%d] traceR", row);
      }
      else EPS_INVALID;
    }

// commands with tag [='...'] or [="..."]
    else if ((c = ungetc(getc(in), in)) == '\'' || c == '"') {
      char c2 = 0;

      if (strcmp(str, "omit") == 0 && (n = fscanf(in, "%*['\"]%" MkString(MAXTAGLEN) "[^'\"]%c", e->tag, &c2)) == 2) {
        e->tag[sizeof e->tag-1] = 0;
        cmd |= eps_omit;
        trace("[%d] omit='%s'", row, e->tag);
        ensure(*e->tag, "invalid empty tag (%s:%d)", option.cfg_file, row);
        ensure(c == c2, "invalid tag quotes (%s:%d)", option.cfg_file, row);
        ensure(!(cmd & (eps_goto | eps_gonum)), "omit tag conflicting with goto (%s:%d)", option.cfg_file, row);
      }
      else if (strcmp(str, "goto") == 0 && (n = fscanf(in, "%*['\"]%" MkString(MAXTAGLEN) "[^'\"]%c", e->tag, &c2)) == 2) {
        e->tag[sizeof e->tag-1] = 0;
        e->num = strtod(e->tag, &end);
        cmd |= !*end ? eps_gonum | eps_istr | eps_nofail : eps_goto;
        trace("[%d] goto='%s'%s", row, e->tag, cmd & eps_gonum ? " (num)" : "");
        ensure(*e->tag, "invalid empty tag (%s:%d)", option.cfg_file, row);
        ensure(c == c2, "invalid tag quotes (%s:%d)", option.cfg_file, row);
        ensure(!(cmd & eps_omit), "goto tag conflicting with omit (%s:%d)", option.cfg_file, row);
      }
      else EPS_INVALID;
    }

// commands with [=...]
    else {
      *buf = 0;
      n = fscanf(in, "%" MkString(MkConcat(1,MAXTAGLEN)) "[^ \t\n\r!#]", buf);
      buf[sizeof buf-1]=0;

      ensure(n == 1 && *buf, "invalid assignment '%s' (%s:%d)", buf, option.cfg_file, row);

      // --- eval registers
      if (*str == 'R' && strchr(buf,'R')) {
        ensure(e->op_n < (short)sizeof e->op, "rule has too many operations (%s:%d)", option.cfg_file, row);

        short dst = strtoul(str+1, &end, 10);
        ensure(dst>=0 && dst < REG_MAX && !*end, "invalid register reference '%s' (%s:%d)", str, option.cfg_file, row);
        ensure(!dst || dst>9, "invalid assignment to read-only register R%d (%s:%d)", dst, option.cfg_file, row);

        char  bop=0;
        short src=0, src2=0;
        bool  pfx = *buf != 'R';
        char  op[2] = { pfx ? *buf : 0, 0 };
        int   n = sscanf(buf+pfx, "R%hd%cR%hd", &src, &bop, &src2);

        if (n == 1) {
          trace("[%d] R%d=%sR%d", row, dst, op, src);
          ensure(!pfx || strchr(REG_UNARY_OP,*buf), "invalid unary operation '%c' (%s:%d)", *buf, option.cfg_file, row);
          e->dst [e->op_n] = dst;
          e->src [e->op_n] = reg_encode(src, *op);
          e->src2[e->op_n] = 0;
          e->op  [e->op_n] = 0;
          e->op_n++;
        }
        else if (n == 3) {
          trace("[%d] R%d=R%d%cR%d", row, dst, src, bop, src2);
          ensure(!pfx, "unallowed prefix operation '%s' (%s:%d)", buf, option.cfg_file, row);
          ensure(strchr(REG_BINARY_OP,bop), "invalid binary operation '%c' (%s:%d)", bop, option.cfg_file, row);
          ensure(reg_isvalid(src ), "invalid register reference '%d' (%s:%d)", src , option.cfg_file, row);
          ensure(reg_isvalid(src2), "invalid register reference '%d' (%s:%d)", src2, option.cfg_file, row);
          e->dst [e->op_n] = dst;
          e->src [e->op_n] = src;
          e->src2[e->op_n] = src2;
          e->op  [e->op_n] = bop;
          e->op_n++;
        }
        else EPS_INVALID;
      }

      // --- load registers
      else if (*str != 'R' && strchr(buf,'R')) {
        bool  pfx = *buf != 'R';
        char  op[2] = { pfx ? *buf : 0, 0 };
        short rn = strtoul(buf+pfx+1, &end, 10);
        ensure(reg_isvalid(rn) && !*end, "invalid register reference '%s' (%s:%d)", buf, option.cfg_file, row);
        ensure(!pfx || strchr(REG_UNARY_OP,*buf), "invalid unary operation '%c' (%s:%d)", *buf, option.cfg_file, row);

             if (strcmp(str, "goto") == 0) {
          cmd |= eps_gonum | eps_istr | eps_nofail;
                            e->gto_reg = reg_encode(rn, *op);  trace("[%d] goto=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "lhs") == 0) {
           cmd |= eps_lhs;  e->lhs_reg = reg_encode(rn, *op);  trace("[%d] lhs=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "rhs") == 0) {
           cmd |= eps_rhs;  e->rhs_reg = reg_encode(rn, *op);  trace("[%d] rhs=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "scl") == 0) {
                            e->scl_reg = reg_encode(rn, *op);  trace("[%d] scl=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "off") == 0) {
                            e->off_reg = reg_encode(rn, *op);  trace("[%d] off=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "abs") == 0) {
          cmd |= eps_abs;   e->abs_reg = reg_encode(rn, *op);  trace("[%d] abs=%sR%d", row, op, rn);  e->_abs_reg = e->abs_reg;
        }
        else if (strcmp(str, "-abs") == 0) {
          cmd |= eps_abs;   e->_abs_reg = reg_encode(rn, *op);  trace("[%d] -abs=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "rel") == 0) {
          cmd |= eps_rel;   e->rel_reg = reg_encode(rn, *op);  trace("[%d] rel=%sR%d", row, op, rn);  e->_rel_reg = e->rel_reg;
        }
        else if (strcmp(str, "-rel") == 0) {
          cmd |= eps_rel;   e->_rel_reg = reg_encode(rn, *op);  trace("[%d] -rel=%sR%d", row, op, rn);
        }
        else if (strcmp(str, "dig") == 0) {
          cmd |= eps_dig;   e->dig_reg = reg_encode(rn, *op);  trace("[%d] dig=%sR%d", row, op, rn);  e->_dig_reg = e->dig_reg;
        }
        else if (strcmp(str, "-dig") == 0) {
          cmd |= eps_dig;   e->_dig_reg = reg_encode(rn, *op);  trace("[%d] -dig=%sR%d", row, op, rn);
        }
        else EPS_INVALID;
      }

      // --- load numbers
      else if (*str != 'R' && !strchr(buf,'R')) {
        double val = strtod(buf, &end);
        ensure(!*end, "invalid number '%s'", buf, option.cfg_file, row);

             if (strcmp(str, "lhs") == 0) {
           cmd |= eps_lhs;  e->lhs = val;  trace("[%d] lhs=%g", row, val);
        }
        else if (strcmp(str, "rhs") == 0) {
           cmd |= eps_rhs;  e->rhs = val;  trace("[%d] rhs=%g", row, val);
        }
        else if (strcmp(str, "scl") == 0) {
                            e->scl = val;  trace("[%d] scl=%g", row, val);
        }
        else if (strcmp(str, "off") == 0) {
                            e->off = val;  trace("[%d] off=%g", row, val);
        }
        else if (strcmp(str, "abs") == 0) {
          cmd |= eps_abs;   e->abs = val;  trace("[%d] abs=%g", row, val);  e->_abs = -val;
          ensure(e->abs >= 0.0 && (cmd & eps_large || e->abs <= 1.0), "invalid absolute constraint (%s:%d)", option.cfg_file, row);
        }
        else if (strcmp(str, "-abs") == 0) {
          cmd |= eps_abs;   e->_abs = val;  trace("[%d] -abs=%g", row, val);
          ensure(e->_abs <= 0.0 && (cmd & eps_large || e->_abs >= -1.0), "invalid negative absolute constraint (%s:%d)", option.cfg_file, row);
        }
        else if (strcmp(str, "rel") == 0) {
          cmd |= eps_rel;   e->rel = val;  trace("[%d] rel=%g", row, val);  e->_rel = -val;
          ensure(e->rel >= 0.0 && (cmd & eps_large || e->rel <= 1.0), "invalid relative constraint (%s:%d)", option.cfg_file, row);
        }
        else if (strcmp(str, "-rel") == 0) {
          cmd |= eps_rel;   e->_rel = val;  trace("[%d] -rel=%g", row, val);
          ensure(e->_rel <= 0.0 && (cmd & eps_large || e->_rel >= -1.0), "invalid negative relative constraint (%s:%d)", option.cfg_file, row);
        }
        else if (strcmp(str, "dig") == 0) {
          cmd |= eps_dig;   e->dig = val;  trace("[%d] dig=%g", row, val);  e->_dig = -val;
          ensure(e->dig >= 1.0, "invalid digital relative constraint (%s:%d)", option.cfg_file, row);
        }
        else if (strcmp(str, "-dig") == 0) {
          cmd |= eps_dig;   e->_dig = val;  trace("[%d] -dig=%g", row, val);
          ensure(e->_dig <= -1.0, "invalid negative digital relative constraint (%s:%d)", option.cfg_file, row);
        }
        else EPS_INVALID;

      }
      else EPS_INVALID;
    }

    // next char
    ungetc((c = getc(in)), in);
    if (c == EOF || (isspace(c) && !isblank(c)) || c == '#' || c == '!') break; 
  }

  // default: abs=eps
  if (!(cmd & (eps_chk | eps_sgg))) {
    cmd |= eps_abs;  e->abs = DBL_MIN;  trace("[%d] abs=%g", row, e->abs);
  }

  // cleanup non-persistant flags (e.g. large)
  e->cmd = (enum eps_cmd)(cmd & eps_mask);  // cast needed because of icc spurious warnings

  trace("<-readEps cmd = %d, str = '%s', c = '%c'", cmd, str, c);

  return cmd == eps_invalid || n == EOF ? EOF : 0;
}

// ----- interface

void
constraint_print(const T* cst, FILE *out)
{
  char op[2] = { 0, 0 };
  short rn=0;

  if (!out) out = stdout;
  if (!cst) { fprintf(out, "(null)"); return; }

  printSlc(&cst->row, out);
  putc(' ', out);
  printSlc(&cst->col, out);
  putc(' ', out);

  if (cst->eps.cmd & eps_alt)    fprintf(out, "alt ");
  if (cst->eps.cmd & eps_any)    fprintf(out, "any ");
  if (cst->eps.cmd & eps_equ)    fprintf(out, "equ ");
  if (cst->eps.cmd & eps_ign)    fprintf(out, "ign ");
  if (cst->eps.cmd & eps_istr)   fprintf(out, "istr ");
  if (cst->eps.cmd & eps_skip)   fprintf(out, "skip ");
  if (cst->eps.cmd & eps_eval)   fprintf(out, "eval ");
  if (cst->eps.cmd & eps_nofail) fprintf(out, "nofail ");
  if (cst->eps.cmd & eps_trace)  fprintf(out, "trace%s ", cst->eps.cmd & eps_traceR ? "R":"");

  if (cst->eps.cmd & eps_omit)   fprintf(out, "omit='%s' ", cst->eps.tag);
  if (cst->eps.cmd & eps_goto)   fprintf(out, "goto='%s' ", cst->eps.tag);

// --- loads
  if (cst->eps.cmd & eps_gonum)  {
    if (cst->eps.gto_reg) {
      rn = reg_decode(cst->eps.gto_reg, op);
      fprintf(out, "goto=%sR%d (num) ", op, rn);
    }
    else fprintf(out, "goto='%s' (num) ", cst->eps.tag);
  }

  if (cst->eps.cmd & eps_lhs) {
    if (cst->eps.lhs_reg) {
      rn = reg_decode(cst->eps.lhs_reg, op);
      fprintf(out, "lhs=%sR%d ", op, rn);
    }
    else fprintf(out, "lhs=%g ", cst->eps.lhs);
  }

  if (cst->eps.cmd & eps_rhs) {
    if (cst->eps.rhs_reg) {
      rn = reg_decode(cst->eps.rhs_reg, op);
      fprintf(out, "rhs=%sR%d ", op, rn);
    }
    else fprintf(out, "rhs=%g ", cst->eps.rhs);
  }

  if (cst->eps.scl_reg) {
    rn = reg_decode(cst->eps.scl_reg, op);
    fprintf(out, "scl=%sR%d ", op, rn);
  }
  else if (cst->eps.scl != 1.0)
    fprintf(out, "scl=%g ", cst->eps.scl);

  if (cst->eps.off_reg) {
    rn = reg_decode(cst->eps.off_reg, op);
    fprintf(out, "off=%sR%d ", op, rn);
  }
  else if (cst->eps.off != 0.0)
    fprintf(out, "off=%g ", cst->eps.off);

  if (cst->eps.cmd & eps_abs) {
    if (cst->eps.abs_reg) {
      rn = reg_decode(cst->eps.abs_reg, op);
      fprintf(out, "abs=%sR%d ", op, rn);
    }
    else {
      if (cst->eps.abs == DBL_MIN)
           fprintf(out, "abs=eps ");
      else fprintf(out, "%sabs=%g ", cst->eps.abs > 1 ? "large " : "", cst->eps.abs);
    }

    if (cst->eps._abs_reg && cst->eps._abs_reg != cst->eps.abs_reg) {
      rn = reg_decode(cst->eps._abs_reg, op);
      fprintf(out, "-abs=%sR%d ", op, rn);
    }
    else if (cst->eps._abs != -cst->eps.abs) {
      if (cst->eps._abs == -DBL_MIN)
           fprintf(out, "-abs=-eps ");
      else fprintf(out, "%s-abs=%g ", cst->eps._abs < -1 ? "large " : "", cst->eps._abs);
    }
  }

  if (cst->eps.cmd & eps_rel) {
    if (cst->eps.rel_reg) {
      rn = reg_decode(cst->eps.rel_reg, op);
      fprintf(out, "rel=%sR%d ", op, rn);
    }
    else fprintf(out, "%srel=%g ", cst->eps.rel > 1 ? "large " : "", cst->eps.rel);

    if (cst->eps._rel_reg && cst->eps._rel_reg != cst->eps.rel_reg) {
      rn = reg_decode(cst->eps._rel_reg, op);
      fprintf(out, "-rel=%sR%d ", op, rn);
    }
    else if (cst->eps._rel != -cst->eps.rel)
      fprintf(out, "%s-rel=%g ", cst->eps._rel < -1 ? "large " : "", cst->eps._rel);
  }

  if (cst->eps.cmd & eps_dig) {
    if (cst->eps.dig_reg) {
      rn = reg_decode(cst->eps.dig_reg, op);
      fprintf(out, "dig=%sR%d ", op, rn);
    }
    else fprintf(out, "dig=%g ", cst->eps.dig);

    if (cst->eps._dig_reg && cst->eps._dig_reg != cst->eps.dig_reg) {
      rn = reg_decode(cst->eps._dig_reg, op);
      fprintf(out, "-dig=%sR%d ", op, rn);
    }
    else if (cst->eps._dig != -cst->eps.dig)
      fprintf(out, "-dig=%g ", cst->eps._dig);
  }   

// --- operations
  for (int j=0; j < cst->eps.op_n; j++) {
    if (cst->eps.op[j])
      fprintf(out, "R%d=R%d%cR%d ", cst->eps.dst[j], cst->eps.src[j], cst->eps.op[j], cst->eps.src2[j]);
    else {
      rn = reg_decode(cst->eps.src[j], op);
      fprintf(out, "R%d=%sR%d ", cst->eps.dst[j], op, rn);
    }
  }
}

void
constraint_scan(T* cst, FILE *in, int *row)
{
  int c;
  assert(cst && row);

  *cst = (T){ .eps = { .cmd = eps_invalid, .scl=1.0 } };

  if (!in) in = stdin;

retry:

  while((c = getc(in)) != EOF && isblank(c)) ;

  // end of file
  if (c == EOF) return;

  ungetc(c, in);

  // comment or empty line
  if (c == '\n' || c == '\r' || c == '#' || c == '!') {
    if (skipLine(in, 0) == '\n') ++*row;
    goto retry;
  }

  cst->idx  = -1;
  cst->line = *row;

  ensure(readSlcOrRng(&cst->row, in      ) != EOF, "invalid row range (%s:%d)"   , option.cfg_file, *row);
  ensure(readSlcOrRng(&cst->col, in      ) != EOF, "invalid column range (%s:%d)", option.cfg_file, *row);
  ensure(readEps     (&cst->eps, in, *row) != EOF, "invalid constraint or command (%s:%d)", option.cfg_file, *row);

  // expand to all columns
  if (cst->eps.cmd & eps_skip || cst->eps.cmd & eps_goto)
    cst->col = slice_initAll();

  // adjust row count
  if (skipLine(in, 0) == '\n') ++*row;
}

