#include <stk_expreval/Variable.hpp>

namespace stk {
namespace expreval {

VariableMap::Resolver &
VariableMap::getDefaultResolver()
{
  static DefaultResolver default_resolver;

  return default_resolver;
}

} // namespace expreval
} // namespace stk
