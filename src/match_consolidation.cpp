#include "match_consolidation.hpp"

// Join adjacent matches that overlap two segments
void consolidate_overlap_match(std::vector<stellar_match> & matches, stellar_match const & match)
{
    auto find_adjacent_match = [&match](stellar_match const & other)
    {
        return match.ref_loc_overlap(other) &&
                match.is_forward_match == other.is_forward_match &&
                match.query_loc_overlap(other);
    };

    auto it = std::find_if(matches.begin(), matches.end(), find_adjacent_match);

    if (it == matches.end())
    {
        matches.push_back(match);
        return;
    }

    while (it != matches.end())
    {
        /* Case A
        segments |------------|
                           |------------|
        match           xxxxxx
        other               xxxxxx

        joined          xxxxxxxxx
        */
        if (match.dbegin < (*it).dbegin)
        {
            stellar_match match_copy = *it;
            // switch this and match and reuse code from Case B
            *it = match;
            (*it).join_adjacent_matches(match_copy);
        }

        /* Case B
        segments |------------|
                           |------------|
        match               xxxxxx
        other           xxxxxx

        joined          xxxxxxxxx
        */

        else if (match.dbegin > (*it).dbegin)
        {
            (*it).join_adjacent_matches(match);
        }

        /* Case C
        segments |------------|
                        |------------|
        match            xxxxx
        other            xxxxx

        */
        else if (match.dend == (*it).dend)
        {
            if (match.percid > (*it).percid)
            {
                *it = match;
            }
        }

        it = std::find_if(++it, matches.end(), find_adjacent_match);
    }
}

void write_output_gff(std::filesystem::path const & out_path, std::vector<stellar_match> const & matches, bool append = false)
{
    std::ofstream fout;

    if (append)
        fout.open(out_path, std::ios_base::app);
    else
        fout.open(out_path);

    for (auto match : matches)
        fout << match.to_string();

    fout.close();
}

void process_matches(consolidation_arguments const & arguments)
{
    std::ifstream fin(arguments.input_file);

    std::string line;
    std::vector<stellar_match> matches;
    std::vector<stellar_match> overlap_matches;
    while (getline(fin, line))
    {
        auto line_vec = get_line_vector<std::string>(line, '\t');
        assert(line_vec.size() == 9); // Stellar GFF format output has 9 columns
        stellar_match match(line_vec, arguments.overlap_length);

        if (match.is_in_segment_overlap())
            consolidate_overlap_match(overlap_matches, match);
        else
            matches.push_back(match);
    }

    fin.close();

    write_output_gff(arguments.output_file, matches);
    write_output_gff(arguments.output_file, overlap_matches, true);
}
