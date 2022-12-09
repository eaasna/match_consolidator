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

struct consolidation_arguments
{
    std::filesystem::path input_file{};
    std::filesystem::path output_file{};
    bool verbose{false};
};

struct stellar_match
{
    std::string dname{};
    size_t segment_start{};
    size_t segment_length{};

    uint64_t dbegin{};
    uint64_t dend{};
    float percid{};
    bool is_forward_match{true};

    std::string query_id{};
    uint64_t qbegin{};
    uint64_t qend{};
    std::string cigar{};
    std::string mutations{};

    // Stellar GFF attributes
    // 1;seq2Range=1280,1378;cigar=97M1D2M;mutations=14A,45G,58T,92C
    std::vector<std::string> attribute_keys{"seq2Range=", "cigar=", "mutations="};

    stellar_match(std::vector<std::string> const match_vec)
    {
        size_t first_delim = match_vec[0].find("_");
        size_t second_delim = match_vec[0].find("_", first_delim + 1);

        dname = match_vec[0].substr(0, first_delim);
        segment_start = stoi(match_vec[0].substr(first_delim + 1, second_delim - first_delim));
        segment_length = stoi(match_vec[0].substr(second_delim + 1));

        dbegin = stoi(match_vec[3]) + segment_start;
        dend = stoi(match_vec[4]) + segment_start;
        percid = std::stof(match_vec[5]);

        if (match_vec[6] == "-")
            is_forward_match = false;

        auto attribute_vec = get_line_vector<std::string>(match_vec[8], ';');
        query_id = attribute_vec[0];
        uint8_t delim_pos = attribute_vec[1].find(",");
        qbegin = stoi(attribute_vec[1].substr(attribute_keys[0].size(), delim_pos - attribute_keys[0].size()));
        qend = stoi(attribute_vec[1].substr(delim_pos + 1));
        cigar = attribute_vec[2].substr(attribute_keys[1].size());
        mutations = attribute_vec[3].substr(attribute_keys[2].size());
    }

    std::string get_percid_str()
    {
        std::stringstream stream;
        // Stellar outputs floating point percentages with 4 decimal precision
        stream << std::fixed << std::setprecision(4) << percid;
        std::string str_percid = stream.str();
        // Remove trailing 0s
        str_percid.erase ( str_percid.find_last_not_of('0') + 1, std::string::npos );
        str_percid.erase ( str_percid.find_last_not_of('.') + 1, std::string::npos );

        return str_percid;
    }

    std::string to_string()
    {
        std::string match_str = dname;
        match_str += "\tStellar\teps-matches\t";
        match_str += std::to_string(dbegin);
        match_str += "\t";
        match_str += std::to_string(dend);
        match_str += "\t";
        match_str += get_percid_str();

        match_str += "\t";

        if (is_forward_match)
            match_str += "+";
        else
            match_str += "-";

        match_str += "\t.\t";
        match_str += query_id;
        match_str += ";";
        match_str += attribute_keys[0];
        match_str += std::to_string(qbegin);
        match_str += ",";
        match_str += std::to_string(qend);
        match_str += ";";
        match_str += attribute_keys[1];
        match_str += cigar;
        match_str += ";";
        match_str += attribute_keys[2];
        match_str += mutations;
        match_str += "\n";

        return match_str;
    }

    bool operator==(stellar_match const & other)
    {
        if (dname == other.dname && dbegin == other.dbegin)
            return true;
        else
            return false;
    }
};