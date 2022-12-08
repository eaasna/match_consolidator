#include "match_consolidation.hpp"

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
            if ((*duplicate_loc).dend < match.dend)
            {
                matches.erase(duplicate_loc);
                matches.push_back(match);
            }
        }
    }

    for (auto & match : matches)
        match.print();

    fin.close();
}
