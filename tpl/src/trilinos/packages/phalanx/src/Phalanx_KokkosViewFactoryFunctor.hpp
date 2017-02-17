#ifndef KOKKOS_VIEW_FACTORY_FUNCTOR_HPP
#define KOKKOS_VIEW_FACTORY_FUNCTOR_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_any.hpp"
#include <string>
#include <typeinfo>

namespace PHX {

  template<typename EvalT>
  class KokkosViewFactoryFunctor {
    
    std::unordered_map<std::string,PHX::any>& fields_;
    const PHX::FieldTag& tag_;
    const std::vector<PHX::index_size_type> extended_dimensions_;
    
  public:
    
    KokkosViewFactoryFunctor(std::unordered_map<std::string,PHX::any>& fields,
			     const PHX::FieldTag& tag,
			     const std::vector<PHX::index_size_type>& extended_dimensions) :
      fields_(fields),
      tag_(tag),
      extended_dimensions_(extended_dimensions)
    {}
    
    template <typename ScalarT>
    void operator()(ScalarT t) const
    {
      if (tag_.dataTypeInfo() == typeid(t)) {
	PHX::KokkosViewFactory<ScalarT,PHX::Device> factory;
	fields_[tag_.identifier()] = factory.buildView(tag_,extended_dimensions_);
      }
    }
    
  };

}

#endif
