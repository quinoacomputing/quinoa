#ifndef MUELU_PARAMETERLISTACCEPTOR_HPP
#define MUELU_PARAMETERLISTACCEPTOR_HPP

#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// TODO See also: Teuchos::ParameterListAcceptor, Teko::Clonable

namespace MueLu {
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::RCP;

  // Printing the valid parameter list only showing documentation fields
  inline void printParameterListOptions(std::ostream &os, const Teuchos::ParameterList & p) {
    p.print(os, Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));
    os << std::endl;
  }

  /*! @class ParameterListAcceptor
      @brief Abstract interface of a class accepting parameter lists.

  */
  class ParameterListAcceptor {

  public:

    //! @name Constructors/Destructors
    //! @{
    ParameterListAcceptor() { }

    virtual ~ParameterListAcceptor() { }

    //! @}
    //! @name Get/Set methods.
    //! @{

    /*! @brief Return a const parameter list of valid parameters that setParameterList() will accept.

       Also define the default values of parameters according to the input parameter list.
    */
    // Questions:
    // - Do we want this function to be static?
    // => yes, but we also want it virtual. So static is not possible. But this method does not access any instance specific data.
    // - Do we want to pass an optional param list as input parameter?
    //   yes => if some parameters depends of other parameters, this would allow to return different default / valid parameters according to already set parameters.
    //   Ex: fact: ILUT => extra param threshold; fact: ILUK => extra param level-of-fill
    //
    // If a parameter is unused AND default, it is printed as [default] by std::cout << paramList but printed as unused by paramList.unused(std::cout).
    // So if we set default parameters in getValidParameters() that are not used, user will get a warning message. We don't want that for [default].
    // One solution is to never set any unused parameter in getValidParameters().
    // If some parameters are available only conditionnaly, do not set them by default in setValidParameters when the conditions are not met.
    //
    virtual RCP<const ParameterList> GetValidParameterList() const = 0;

    //! @brief Set parameters from a parameter list and return with default values.
    //
    // Questions:
    //  - Input: const ParameterList or not const? IE: do we want to modify user paramlist directly?
    //    => not const: user can make a copy outside if he do not want us to modify his initial paramlist.
    //  - Input: RCP or just reference? RCP avoid the copy of the parameter list but
    //           I'm afraid that a user modify the list outside of the class after doing a SetParameterList
    //           (this by-pass the validation and can create some confusion).
    //           In another hand, if we use a copy of the list, we do not set the "[used]"/["unused"] flag
    //           on the original list during Build(). GetParameterList has to be called to retrieve the list with the correct flags.
    //
    // What we really want is for the user to have a const version outside and MueLu having a non-const version inside.
    // Is there a C++ pattern to do that?
    virtual void SetParameterList(const ParameterList& paramList) = 0;

    // Get the parameter list.
    // Teuchos::ParameterListAcceptor also provides a non-const version of this but I don't see why.
    // We don't want a user to mess with the internal parameter list.
    // To change a parameter, one can make a copy and call SetParameterList again.
    // No need for RCP, it's just a view
    virtual const Teuchos::ParameterList & GetParameterList() const = 0;

    //! @brief Set a parameter directly as a ParameterEntry.
    virtual void SetParameter(const std::string &name, const ParameterEntry &entry) = 0;

    //! @brief Retrieves a const entry with the name name.
    virtual const ParameterEntry & GetParameter(const std::string &name) const = 0;

    virtual void GetDocumentation(std::ostream &os) const = 0;

    //! @}

  }; //class ParameterListAcceptor

  // Partial implementation of ParameterListAcceptor that stores the object parameters in an internal parameterList
  class ParameterListAcceptorImpl: public ParameterListAcceptor {

  public:

    ParameterListAcceptorImpl() { }

    virtual ~ParameterListAcceptorImpl() {
      // TAW: 3/17/2016
      // The following code is useless. Comment it out to make coverity happy (CID22611)
      /*bool warnings = false; //TODO
      if (warnings)
        paramList_.unused(std::cout);*/
    }

    virtual void SetParameterList(const ParameterList& paramList) {
      // This call is only for cosmetic reasons.
      // If one calls SetParameterList before GetParameterList, that would mean
      // that paramList_ has not been initialized yet. Therefore, the parameter
      // would be put in it in the order user provided, and not in the order of
      // valid parameter list. We'd like to have consistency in the order, so
      // we do this extra call, which is no-op if paramList_ has already been
      // initialized.
      paramList_ = GetParameterList();

      paramList_.setParameters(paramList);

      // Validate and add defaults parameters.
      RCP<const ParameterList> validParamList = GetValidParameterList();
      if (validParamList != Teuchos::null) {
        paramList_.validateParametersAndSetDefaults(*validParamList);
      } else {
        // Teuchos::null means that GetValidParameterList() not implemented.
        // As such, we skip validation and have not way to set default values,
        // which is potentially dangerous
      }
    }

    // The returned list always has an entry for each valid parameter.
    // Therefore, there is not need to test if a parameter is present before getting it.
    virtual const Teuchos::ParameterList& GetParameterList() const {
      if (paramList_.numParams() == 0) {
        // Set paramList_ to the default list
        RCP<const ParameterList> validParamList = GetValidParameterList();
        if (validParamList != Teuchos::null) {
          // Instead of simply doing
          //   paramList_ = *validParamList;
          // we use more complicated Teuchos calls, because we would like to
          // have [default] values in the beginning
          paramList_.validateParametersAndSetDefaults(*validParamList);
        }

      } else {
        // We are sure that the list has all the valid parameters defined
        // because the parameter list validation process adds the default
        // values to the user list
      }

      return paramList_;
    }

    void SetParameter(const std::string& name, const ParameterEntry& entry) {
      Teuchos::ParameterList paramList;
      paramList.setEntry(name, entry);
      SetParameterList(paramList); // This forces revalidation of the list
    }

    const ParameterEntry & GetParameter(const std::string &name) const {
      return GetParameterList().getEntry(name);
    }

    virtual void GetDocumentation(std::ostream &os) const {
      // default implementation

      RCP<const ParameterList> validParamList = GetValidParameterList();
      if (validParamList == Teuchos::null) {
        os << "## Documentation not available:" << std::endl;
        return;
      }

      os << "## Parameters:" << std::endl;
      printParameterListOptions(os, *validParamList);

      os << "## Fully described default method:" << std::endl;
      validParamList->print(os, 2, true, false);
      os << std::endl;
    }

  private:

    mutable // mutable can be avoid by changing return type of GetParameterList() to RCP but conceptually, I like the fact that GetParameterList returns ref to param list.
    Teuchos::ParameterList paramList_; // Private: Use GetParameterList() to access this list
                                       // This list might be empty before calling GetParameterList()
                                       // but it is populate with automatically if needed in SetParameterList/GetParameterList/...
                                       // Therefore, there is no need to test if a parameter exist before accessing it with pL.get("paramName")
                                       // in sub classes.

    // Note: Teuchos::ParameterListAcceptor has a getMyParamList() function to access paramList_ without going through the logic of getParameterList()
    // Do we need that too?

  };


}

#endif // MUELU_PARAMETERLISTACCEPTOR_HPP
