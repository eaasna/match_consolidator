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

struct gff_mutations
{
    gff_mutations() = default;
    gff_mutations(gff_mutations const &) = default;
    gff_mutations & operator=(gff_mutations const &) = default;
    gff_mutations(gff_mutations &&) = default;
    gff_mutations & operator=(gff_mutations &&) = default;
    ~gff_mutations() = default;

    struct mutation
    {
        uint16_t pos;
        char allele;

        mutation(std::string mutation_str)
        {
            pos = std::stoi(mutation_str.substr(0, mutation_str.size() - 1));
            allele = mutation_str[mutation_str.size() - 1];
        }

        void increase_pos(uint16_t const increase)
        {
            pos += increase;
        }
    };

    uint64_t abs_pos;
    uint16_t max_pos{0};
    std::vector<mutation> mut_vec;
    std::string prefix = "mutations=";

    gff_mutations(std::string const mutations_str, uint64_t const dbegin)
    {
        std::istringstream iss(mutations_str.substr(prefix.size()));
        std::string field;
        while(std::getline(iss, field, ','))
        {
            mutation mut(field);
            mut_vec.emplace_back(mut);
            if (mut.pos > max_pos)
                max_pos = mut.pos;
        }

        abs_pos = dbegin;
    }

    void join_mutations(gff_mutations const & other)
    {

        //!TODO: need to add this* mutations to the end of other* mutations
        auto mut_it = std::find_if(other.mut_vec.begin(),
                                other.mut_vec.end(),
                                [&](const auto & mut) { return mut.pos > other.abs_pos + max_pos; });

        for (;mut_it != other.mut_vec.end(); mut_it++)
        {
            mutation adjusted_mut = *mut_it;
            adjusted_mut.increase_pos(other.abs_pos - abs_pos);
            mut_vec.push_back(adjusted_mut);
        }
    }

    std::string to_string()
    {
        std::string mut_str = prefix;

        for (auto mut : mut_vec)
        {
            mut_str += std::to_string(mut.pos);
            mut_str += mut.allele;
            mut_str += ",";
        }

        return mut_str.substr(0, mut_str.size() - 1);
    }

};

struct stellar_match
{
    std::string dname{};
    size_t segment_start{};
    size_t segment_length{};
    size_t segment_overlap{};

    uint64_t dbegin{};
    uint64_t dend{};
    float percid{};
    bool is_forward_match{true};

    std::string query_id{};
    uint64_t qbegin{};
    uint64_t qend{};
    std::string cigar{};
    gff_mutations mutations;

    // Stellar GFF attributes
    // 1;seq2Range=1280,1378;cigar=97M1D2M;mutations=14A,45G,58T,92C
    std::vector<std::string> attribute_keys{"seq2Range=", "cigar=", "mutations="};

    stellar_match(std::vector<std::string> const match_vec, uint16_t const seg_overlap)
    {
        size_t first_delim = match_vec[0].find("_");
        size_t second_delim = match_vec[0].find("_", first_delim + 1);
        segment_overlap = seg_overlap;

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

        mutations = gff_mutations(attribute_vec[3], dbegin);
    }

    bool is_in_segment_overlap()
    {
        if (dbegin <= segment_start + segment_overlap ||
            dend >= segment_start + segment_length - segment_overlap)
            return true;

        return false;
    }

    bool operator==(stellar_match const & other)
    {
        if (dname == other.dname && (dbegin == other.dbegin || dend == other.dend))
            return true;
        else
            return false;
    }

    bool ref_loc_overlap(stellar_match const & other) const
    {
        if (dname == other.dname)
        {
            if (ranges_overlap(dbegin, dend, other.dbegin, other.dend))
                return true;
        }
        return false;
    }

    bool query_loc_overlap(stellar_match const & other) const
    {
        if (query_id == other.query_id)
        {
            if (ranges_overlap(qbegin, qend, other.qbegin, other.qend))
                return true;
        }
        return false;
    }

    // TODO: update percid, cigar
    /* Case A
    segments |------------|
                       |------------|
    match           xxxxxx
    this                xxxxxx

    joined          xxxxxxxxx
    */
    void join_adjacent_matches(stellar_match const & match)
    {
        dbegin = match.dbegin;
        qbegin = match.qbegin;
        percid = 100;

        mutations.join_mutations(match.mutations);
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
        match_str += mutations.to_string();
        match_str += "\n";

        return match_str;
    }

};