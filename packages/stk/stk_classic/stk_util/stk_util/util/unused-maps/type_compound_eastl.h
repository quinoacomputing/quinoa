/*
Copyright (C) 2005,2009-2010 Electronic Arts, Inc.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1.  Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
3.  Neither the name of Electronic Arts, Inc. ("EA") nor the names of
    its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY ELECTRONIC ARTS AND ITS CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ELECTRONIC ARTS OR ITS CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

///////////////////////////////////////////////////////////////////////////////
// EASTL/internal/type_compound.h
// Written and maintained by Paul Pedriana - 2005.
///////////////////////////////////////////////////////////////////////////////


#ifndef EASTL_INTERNAL_TYPE_COMPOUND_H
#define EASTL_INTERNAL_TYPE_COMPOUND_H


namespace eastl
{

    // The following properties or relations are defined here. If the given
    // item is missing then it simply hasn't been implemented, at least not yet.
    //   is_array
    //   is_pointer
    //   is_reference
    //   is_member_object_pointer
    //   is_member_function_pointer
    //   is_member_pointer
    //   is_enum
    //   is_union
    //   is_class
    //   is_polymorphic
    //   is_function
    //   is_object
    //   is_scalar
    //   is_compound
    //   is_same
    //   is_convertible

    ///////////////////////////////////////////////////////////////////////
    // is_array
    //
    // is_array<T>::value == true if and only if T is an array type.
    // As of this writing, the SNC compiler (EDG-based) doesn't compile
    // the code below and says that returning an array is illegal.
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    T (*is_array_tester1(empty<T>))(empty<T>);
    char is_array_tester1(...);     // May need to use __cdecl under VC++.

    template <typename T>
    no_type  is_array_tester2(T(*)(empty<T>));
    yes_type is_array_tester2(...); // May need to use __cdecl under VC++.

    template <typename T>
    struct is_array_helper {
        static empty<T> emptyInstance;
    };

    template <typename T>
    struct is_array : public integral_constant<bool,
        sizeof(is_array_tester2(is_array_tester1(is_array_helper<T>::emptyInstance))) == 1
    >{};



    ///////////////////////////////////////////////////////////////////////
    // is_reference
    //
    // is_reference<T>::value == true if and only if T is a reference type.
    // This category includes reference to function types.
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_reference : public false_type{};
    template <typename T> struct is_reference<T&> : public true_type{};



    ///////////////////////////////////////////////////////////////////////
    // is_member_function_pointer
    //
    // is_member_function_pointer<T>::value == true if and only if T is a
    // pointer to member function type.
    //
    ///////////////////////////////////////////////////////////////////////
    // We detect member functions with 0 to N arguments. We can extend this
    // for additional arguments if necessary.
    // To do: Make volatile and const volatile versions of these in addition to non-const and const.
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_mem_fun_pointer_value : public false_type{};
    template <typename R, typename T> struct is_mem_fun_pointer_value<R (T::*)()> : public true_type{};
    template <typename R, typename T> struct is_mem_fun_pointer_value<R (T::*)() const> : public true_type{};
    template <typename R, typename T, typename Arg0> struct is_mem_fun_pointer_value<R (T::*)(Arg0)> : public true_type{};
    template <typename R, typename T, typename Arg0> struct is_mem_fun_pointer_value<R (T::*)(Arg0) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6) const> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)> : public true_type{};
    template <typename R, typename T, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7> struct is_mem_fun_pointer_value<R (T::*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7) const> : public true_type{};

    template <typename T>
    struct is_member_function_pointer : public integral_constant<bool, is_mem_fun_pointer_value<T>::value>{};


    ///////////////////////////////////////////////////////////////////////
    // is_member_pointer
    //
    // is_member_pointer<T>::value == true if and only if:
    //    is_member_object_pointer<T>::value == true, or
    //    is_member_function_pointer<T>::value == true
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_member_pointer : public integral_constant<bool, is_member_function_pointer<T>::value>{};

    template <typename T, typename U> struct is_member_pointer<U T::*> : public true_type{};




    ///////////////////////////////////////////////////////////////////////
    // is_pointer
    //
    // is_pointer<T>::value == true if and only if T is a pointer type.
    // This category includes function pointer types, but not pointer to
    // member types.
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_pointer_helper : public false_type{};

    template <typename T> struct is_pointer_helper<T*>                : public true_type{};
    template <typename T> struct is_pointer_helper<T* const>          : public true_type{};
    template <typename T> struct is_pointer_helper<T* volatile>       : public true_type{};
    template <typename T> struct is_pointer_helper<T* const volatile> : public true_type{};

    template <typename T>
    struct is_pointer_value : public type_and<is_pointer_helper<T>::value, type_not<is_member_pointer<T>::value>::value> {};

    template <typename T>
    struct is_pointer : public integral_constant<bool, is_pointer_value<T>::value>{};



    ///////////////////////////////////////////////////////////////////////
    // is_same
    //
    // Given two (possibly identical) types T and U, is_same<T, U>::value == true
    // if and only if T and U are the same type.
    //
    ///////////////////////////////////////////////////////////////////////
    template<typename T, typename U>
    struct is_same : public false_type { };

    template<typename T>
    struct is_same<T, T> : public true_type { };


    ///////////////////////////////////////////////////////////////////////
    // is_convertible
    //
    // Given two (possible identical) types From and To, is_convertible<From, To>::value == true
    // if and only if an lvalue of type From can be implicitly converted to type To,
    // or is_void<To>::value == true
    //
    // is_convertible may only be applied to complete types.
    // Type To may not be an abstract type.
    // If the conversion is ambiguous, the program is ill-formed.
    // If either or both of From and To are class types, and the conversion would invoke
    // non-public member functions of either From or To (such as a private constructor of To,
    // or a private conversion operator of From), the program is ill-formed.
    //
    // Note that without compiler help, both is_convertible and is_base
    // can produce compiler errors if the conversion is ambiguous.
    // Example:
    //    struct A {};
    //    struct B : A {};
    //    struct C : A {};
    //    struct D : B, C {};
    //    is_convertible<D*, A*>::value; // Generates compiler error.
    ///////////////////////////////////////////////////////////////////////
    #if !defined(__GNUC__) || (__GNUC__ >= 3) // GCC 2.x doesn't like the code below.
        template <typename From, typename To, bool is_from_void = false, bool is_to_void = false>
        struct is_convertible_helper {
            static yes_type Test(To);  // May need to use __cdecl under VC++.
            static no_type  Test(...); // May need to use __cdecl under VC++.
            static From from;
            typedef integral_constant<bool, sizeof(Test(from)) == sizeof(yes_type)> result;
        };

        // void is not convertible to non-void
        template <typename From, typename To>
        struct is_convertible_helper<From, To, true, false> { typedef false_type result; };

        // Anything is convertible to void
        template <typename From, typename To, bool is_from_void>
        struct is_convertible_helper<From, To, is_from_void, true> { typedef true_type result; };

        template <typename From, typename To>
        struct is_convertible : public is_convertible_helper<From, To, is_void<From>::value, is_void<To>::value>::result {};

    #else
        template <typename From, typename To>
        struct is_convertible : public false_type{};
    #endif


    ///////////////////////////////////////////////////////////////////////
    // is_union
    //
    // is_union<T>::value == true if and only if T is a union type.
    //
    // There is no way to tell if a type is a union without compiler help.
    // As of this writing, only Metrowerks v8+ supports such functionality
    // via 'msl::is_union<T>::value'. The user can force something to be
    // evaluated as a union via EASTL_DECLARE_UNION.
    ///////////////////////////////////////////////////////////////////////
    template <typename T> struct is_union : public false_type{};

    #define EASTL_DECLARE_UNION(T) namespace eastl{ template <> struct is_union<T> : public true_type{}; template <> struct is_union<const T> : public true_type{}; }




    ///////////////////////////////////////////////////////////////////////
    // is_class
    //
    // is_class<T>::value == true if and only if T is a class or struct
    // type (and not a union type).
    //
    // Without specific compiler help, it is not possible to
    // distinguish between unions and classes. As a result, is_class
    // will erroneously evaluate to true for union types.
    ///////////////////////////////////////////////////////////////////////
    #if defined(__MWERKS__)
        // To do: Switch this to use msl_utility type traits.
        template <typename T>
        struct is_class : public false_type{};
    #elif !defined(__GNUC__) || (((__GNUC__ * 100) + __GNUC_MINOR__) >= 304) // Not GCC or GCC 3.4+
        template <typename U> static yes_type is_class_helper(void (U::*)());
        template <typename U> static no_type  is_class_helper(...);

        template <typename T>
        struct is_class : public integral_constant<bool,
            sizeof(is_class_helper<T>(0)) == sizeof(yes_type) && !is_union<T>::value
        >{};
    #else
        // GCC 2.x version, due to GCC being broken.
        template <typename T>
        struct is_class : public false_type{};
    #endif



    ///////////////////////////////////////////////////////////////////////
    // is_enum
    //
    // is_enum<T>::value == true if and only if T is an enumeration type.
    //
    ///////////////////////////////////////////////////////////////////////
    struct int_convertible{ int_convertible(int); };

    template <bool is_arithmetic_or_reference>
    struct is_enum_helper { template <typename T> struct nest : public is_convertible<T, int_convertible>{}; };

    template <>
    struct is_enum_helper<true> { template <typename T> struct nest : public false_type {}; };

    template <typename T>
    struct is_enum_helper2
    {
        typedef type_or<is_arithmetic<T>::value, is_reference<T>::value, is_class<T>::value> selector;
        typedef is_enum_helper<selector::value> helper_t;
        typedef typename add_reference<T>::type ref_t;
        typedef typename helper_t::template nest<ref_t> result;
    };

    template <typename T>
    struct is_enum : public integral_constant<bool, is_enum_helper2<T>::result::value>{};

    template <> struct is_enum<void> : public false_type {};
    template <> struct is_enum<void const> : public false_type {};
    template <> struct is_enum<void volatile> : public false_type {};
    template <> struct is_enum<void const volatile> : public false_type {};

    #define EASTL_DECLARE_ENUM(T) namespace eastl{ template <> struct is_enum<T> : public true_type{}; template <> struct is_enum<const T> : public true_type{}; }


    ///////////////////////////////////////////////////////////////////////
    // is_polymorphic
    //
    // is_polymorphic<T>::value == true if and only if T is a class or struct
    // that declares or inherits a virtual function. is_polymorphic may only
    // be applied to complete types.
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_polymorphic_imp1
    {
        typedef typename remove_cv<T>::type t;

        struct helper_1 : public t
        {
            helper_1();
            ~helper_1() throw();
            char pad[64];
        };

        struct helper_2 : public t
        {
            helper_2();
            virtual ~helper_2() throw();
            #ifndef _MSC_VER
                virtual void foo();
            #endif
            char pad[64];
        };

        static const bool value = (sizeof(helper_1) == sizeof(helper_2));
    };

    template <typename T>
    struct is_polymorphic_imp2{ static const bool value = false; };

    template <bool is_class>
    struct is_polymorphic_selector{ template <typename T> struct rebind{ typedef is_polymorphic_imp2<T> type; }; };

    template <>
    struct is_polymorphic_selector<true>{ template <typename T> struct rebind{ typedef is_polymorphic_imp1<T> type; }; };

    template <typename T>
    struct is_polymorphic_value{
        typedef is_polymorphic_selector<is_class<T>::value> selector;
        typedef typename selector::template rebind<T> binder;
        typedef typename binder::type imp_type;
        static const bool value = imp_type::value;
    };

    template <typename T>
    struct is_polymorphic : public integral_constant<bool, is_polymorphic_value<T>::value>{};




    ///////////////////////////////////////////////////////////////////////
    // is_function
    //
    // is_function<T>::value == true  if and only if T is a function type.
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename R> struct is_function_ptr_helper : public false_type{};
    template <typename R> struct is_function_ptr_helper<R (*)()> : public true_type{};
    template <typename R, typename Arg0> struct is_function_ptr_helper<R (*)(Arg0)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1> struct is_function_ptr_helper<R (*)(Arg0, Arg1)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2, typename Arg3> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2, Arg3)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2, Arg3, Arg4)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6)> : public true_type{};
    template <typename R, typename Arg0, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename Arg6, typename Arg7> struct is_function_ptr_helper<R (*)(Arg0, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7)> : public true_type{};

    template <bool is_ref = true>
    struct is_function_chooser{ template <typename T> struct result_ : public false_type{}; };

    template <>
    struct is_function_chooser<false>{ template <typename T> struct result_ : public is_function_ptr_helper<T*>{}; };

    template <typename T>
    struct is_function_value : public is_function_chooser<is_reference<T>::value>::template result_<T>{};

    template <typename T>
    struct is_function : public integral_constant<bool, is_function_value<T>::value>{};




    ///////////////////////////////////////////////////////////////////////
    // is_object
    //
    // is_object<T>::value == true if and only if:
    //    is_reference<T>::value == false, and
    //    is_function<T>::value == false, and
    //    is_void<T>::value == false
    //
    // The C++ standard, section 3.9p9, states: "An object type is a
    // (possibly cv-qualified) type that is not a function type, not a
    // reference type, and not incomplete (except for an incompletely
    // defined object type).
    ///////////////////////////////////////////////////////////////////////

    template <typename T>
    struct is_object : public integral_constant<bool,
        !is_reference<T>::value && !is_void<T>::value && !is_function<T>::value
    >{};



    ///////////////////////////////////////////////////////////////////////
    // is_scalar
    //
    // is_scalar<T>::value == true if and only if:
    //    is_arithmetic<T>::value == true, or
    //    is_enum<T>::value == true, or
    //    is_pointer<T>::value == true, or
    //    is_member_pointer<T>::value
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_scalar : public integral_constant<bool, is_arithmetic<T>::value || is_enum<T>::value>{};

    template <typename T> struct is_scalar<T*> : public true_type {};
    template <typename T> struct is_scalar<T* const> : public true_type {};
    template <typename T> struct is_scalar<T* volatile> : public true_type {};
    template <typename T> struct is_scalar<T* const volatile> : public true_type {};



    ///////////////////////////////////////////////////////////////////////
    // is_compound
    //
    // Compound means anything but fundamental. See C++ standard, section 3.9.2.
    //
    // is_compound<T>::value == true if and only if:
    //    is_fundamental<T>::value == false
    //
    // Thus, is_compound<T>::value == true if and only if:
    //    is_floating_point<T>::value == false, and
    //    is_integral<T>::value == false, and
    //    is_void<T>::value == false
    //
    ///////////////////////////////////////////////////////////////////////
    template <typename T>
    struct is_compound : public integral_constant<bool, !is_fundamental<T>::value>{};


} // namespace eastl


#endif // Header include guard
