//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Thu Nov 14 08:17:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <list>

namespace tk {

//! Driver base class
class Driver {

  public:
    //! Constructor
    explicit Driver() = default;

    //! Destructor
    virtual ~Driver() noexcept = default;

    //! Execute
    virtual void execute() const = 0;

    //! Register into factory
    //! \param[in] C       Type of the (derived) class constructor
    //! \param[in] F       Type of factory to add to
    //! \param[in] O       Type of option to add
    //! \param[in] E       Type of enum to add
    //! \param[in] Args... Types of variable number of arguments to constructor
    //! \param[in] f       Factory instance to add to
    //! \param[in] reg     List of enums to add enum to
    //! \param[in] o       Type of option to add
    //! \param[in] e       Enum key to factory's std::map
    //! \param[in] args    Variable number of arguments to constructor
    template< class C, class F, class O, typename E, typename... Args >
    void add( F& f, std::list<E>& reg, const O& o, E e, const Args&... args ) {
      reg.push_back( o.template add<C>( f, e, std::move(args)... ) );
    }

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;
};

} // namespace tk

#endif // Driver_h
