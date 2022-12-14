#include "match_consolidation.hpp"

void write_output_gff(std::filesystem::path const & out_path, std::vector<stellar_match> const & matches)
{
    std::ofstream fout(out_path);

    if(fout.is_open())
    {
        for (auto match : matches)
            fout << match.to_string();

        fout.close();
    }
}

void consolidate_matches(consolidation_arguments const & arguments)
{
    std::ifstream fin(arguments.input_file);

    std::string line;
    std::vector<stellar_match> matches;
    while (getline(fin, line))
    {
        // seqan3::debug_stream << line << "\n";
        auto line_vec = get_line_vector<std::string>(line, '\t');
        assert(line_vec.size() == 9); // Stellar GFF format output has 9 columns
        stellar_match match(line_vec, arguments.overlap_length);

        /*
        segments |------------|
                          |------------|
        match           xxxxxx
        other              xxxxxx

        joined          xxxxxxxxx
        */
        auto find_next_segment_match = [&match](const stellar_match & other)
        {
            return match.on_previous_segment(other) &&
                   match.dend >= other.dbegin &&
                   match.is_forward_match == other.is_forward_match &&
                   match.overlaps_query(other);
        };

        for (auto it = std::find_if(matches.begin(), matches.end(), find_next_segment_match);
                  it != matches.end();
                  it = std::find_if(++it, matches.end(), find_next_segment_match))
        {
            (*it).dbegin = match.dbegin;
            (*it).qbegin = match.qbegin;
            (*it).percid = 100;
            // TODO: update percid, cigar, mutations
        }

        /*
        segments |------------|
                          |------------|
        match              xxxxxx
        other           xxxxxx

        joined          xxxxxxxxx
        */
        auto find_previous_segment_match = [&match](const stellar_match & other)
        {
            return match.on_next_segment(other) &&
                   other.dend >= match.dbegin &&
                   match.is_forward_match == other.is_forward_match &&
                   match.overlaps_query(other);
        };

        for (auto it = std::find_if(matches.begin(), matches.end(), find_previous_segment_match);
                  it != matches.end();
                  it = std::find_if(++it, matches.end(), find_previous_segment_match))
        {
            (*it).dend = match.dend;
            (*it).qend = match.qend;
            (*it).percid = 100;
            // TODO: update percid, cigar, mutations
        }

        // remove duplicates on segment overlaps
        auto duplicate_loc = std::find(matches.begin(), matches.end(), match);
        if (duplicate_loc == matches.end())
        {
            matches.push_back(match);
        }
        else
        {
            if ((*duplicate_loc).percid < match.percid)
            {
                matches.erase(duplicate_loc);
                matches.push_back(match);
            }
        }

    }

    fin.close();

    write_output_gff(arguments.output_file, matches);
}
