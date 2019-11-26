#pragma once

#include <iostream>

void print()
{
    std::cout << std::endl;
}

template <typename type>
void print(const type& arg)
{
    std::cout << arg << std::endl;
}

template <typename first, typename... other>
void print(const first& first_arg, const other&... other_args)
{
    std::cout << first_arg << " ";
    print(other_args...);
}