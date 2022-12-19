#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>

template <typename field_t>
std::vector<field_t> get_line_vector(std::string const line, char const delim)
{
    std::vector<field_t> line_vec;

    std::istringstream iss(line);
    std::string field;
    while(std::getline(iss, field, delim))
        line_vec.push_back(field);

    return line_vec;
}

// test whether ranges [a1;a2] and [b1;b2] overlap
template <typename field_t>
inline bool ranges_overlap(field_t a1, field_t a2, field_t b1, field_t b2)
{
    if (a1 <= b2 && b1 <= a2)
        return true;

    return false;
}

struct consolidation_arguments
{
    std::filesystem::path input_file{};
    std::filesystem::path output_file{};
    size_t overlap_length{};
    bool verbose{false};
};
