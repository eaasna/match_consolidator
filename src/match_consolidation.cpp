#include "match_consolidation.hpp"

void consolidate_matches(consolidation_arguments const & arguments)
{
    std::ifstream fin(arguments.input_file);

    std::string line;
    std::vector<stellar_match> matches;
    while (getline(fin, line))
    {
        seqan3::debug_stream << line << "\n";
        auto match_vec = get_line_vector<std::string>(line, '\t');
        assert(match_vec.size() == 9); // Stellar GFF format output has 9 columns
        stellar_match match(match_vec);
        matches.push_back(match);
        match.print();
    }

    fin.close();
}
