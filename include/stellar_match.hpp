#pragma once

#include "gff_mutations.hpp"
#include "gff_percid.hpp"

struct stellar_match
{
    std::string dname{};
    size_t segment_start{};
    size_t segment_length{};
    size_t segment_overlap{};

    uint64_t dbegin{};
    uint64_t dend{};
    gff_percid percid;
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

        percid = gff_percid(match_vec[5]);

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

    // TODO: update cigar
    /* Case A
    segments |------------|
                       |------------|
    match           xxxxxx
    this                xxxxxx

    joined          xxxxxxxxxx
    */
    void join_adjacent_matches(stellar_match const & match)
    {
        dend = match.dend;
        qend = match.qend;
        mutations.join_mutations(match.mutations);
        percid.update(dbegin, mutations, dend);
    }

    std::string to_string()
    {
        std::string match_str = dname;
        match_str += "\tStellar\teps-matches\t";
        match_str += std::to_string(dbegin);
        match_str += "\t";
        match_str += std::to_string(dend);
        match_str += "\t";
        match_str += percid.to_string();

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
