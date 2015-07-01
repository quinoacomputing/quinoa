// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Trilinos includes
#include "Teuchos_ParameterEntry.hpp" // class data element 
#include "Teuchos_Assert.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_map.hpp"
#include "PyTeuchos_Utils.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

Teuchos::ParameterList dict2ParameterList(PyObject* obj);

namespace Teuchos
{
  class ParameterList
  {
  public:
    ParameterList();
    ~ParameterList();

    %extend
    {
      ParameterList(PyObject* dict)
        {
          ParameterList* p = new ParameterList(dict2ParameterList(dict));
          return p;
        }
    }


    bool isParameter(const std::string& name) const ;


    /* setting parameters */

    %extend
    {
      void setString(const std::string& name, const std::string& val)
      {
        self->set<std::string>(name, val);
      }
    }

    %extend
    {
      void setDouble(const std::string& name, const double& val)
      {
        self->set<double>(name, val);
      }
    }

    %extend
    {
      void setInt(const std::string& name, const int& val)
      {
        self->set<int>(name, val);
      }
    }

    %extend
    {
      void setSublist(const std::string& name, const ParameterList& val)
      {
        self->set<Teuchos::ParameterList>(name, val);
      }
    }

    /* checking parameter types */

    %extend
    {
      bool isString(const std::string& name) const 
      {
        return self->isType<std::string>(name);
      }
    }

    %extend
    {
      bool isTypeDouble(const std::string& name) const 
      {
        return self->isType<double>(name);
      }
    }

    %extend
    {
      bool isTypeInt(const std::string& name) const 
      {
        return self->isType<int>(name);
      }
    }

    %extend
    {
      bool isSublist(const std::string& name) const 
      {
        return self->isType<Teuchos::ParameterList>(name);
      }
    }



    /* checking parameter types */

    %extend
    {
      bool isString(const std::string& name) const 
      {
        return self->isType<std::string>(name);
      }
    }

    %extend
    {
      bool isDouble(const std::string& name) const 
      {
        return self->isType<double>(name);
      }
    }

    %extend
    {
      bool isInt(const std::string& name) const 
      {
        return self->isType<int>(name);
      }
    }

    %extend
    {
      bool isSublist(const std::string& name) const 
      {
        return self->isType<ParameterList>(name);
      }
    }


    /* getting parameters */

    %extend
    {
      std::string getString(const std::string& name) const 
      {
        return self->get<std::string>(name);
      }
    }

    %extend
    {
      double getDouble(const std::string& name) const 
      {
        return self->get<double>(name);
      }
    }

    %extend
    {
      int getInt(const std::string& name) const 
      {
        return self->get<int>(name);
      }
    }

    %extend
    {
      const Teuchos::ParameterList getSublist(const std::string& name) const 
      {
        return self->get<Teuchos::ParameterList>(name);
      }
    }



  };

  %extend ParameterList
  {
    using namespace std;
    std::string __str__() 
    {
      std::string rtn; 
      std::stringstream os;
      self->print(os);
      rtn = os.str();
      return rtn;
    }
  }

}





