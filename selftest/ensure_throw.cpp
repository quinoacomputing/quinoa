#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>
#include <string>
#include <stdexcept>

namespace tut
{

/**
 * Testing skip() method.
 */
struct throw_test
{
    virtual ~throw_test()
    {
    }
};

typedef test_group<throw_test> tf;
typedef tf::object object;
tf throw_test("throw");

namespace 
{
    struct foo_exception {};
    struct bar_exception {};

    void foo()
    {
        throw foo_exception();
    }
}

template<>
template<>
void object::test<1>()
{
    set_test_name("checks throw");

    try
    {
        ensure_THROW( foo(), bar_exception );
        throw std::runtime_error("throw doesn't work");
    }
    catch (const failure& ex)
    {
        if (std::string(ex.what()).find("foo()") == std::string::npos )
        {
            fail("throw doesn't contain proper message");
        }
    }
}

template<>
template<>
void object::test<2>()
{
    set_test_name("checks throw");

    try
    {
        ensure_THROW( foo(), foo_exception );
    }
    catch (const failure& ex)
    {
        fail("positive throw doesn't work");
    }
}

template<>
template<>
void object::test<3>()
{
    set_test_name("checks no_throw");

    try
    {
        ensure_NO_THROW( skip() );
        throw std::runtime_error("no_throw doesn't work");
    }
    catch (const failure& ex)
    {
        if (std::string(ex.what()).find("skip()") == std::string::npos )
        {
            fail("no_throw doesn't contain proper message");
        }
    }
}

}

