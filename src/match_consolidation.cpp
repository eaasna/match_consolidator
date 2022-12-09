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
        stellar_match match(line_vec);

        auto duplicate_loc = std::find(matches.begin(), matches.end(), match);
        if (duplicate_loc == matches.end())
        {
            matches.push_back(match);
        }
        else
        {
            // TODO: unify adjacent matches if they overlap two segments
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
