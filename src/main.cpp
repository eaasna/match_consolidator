#include <sstream>

#include <seqan3/argument_parser/all.hpp>

#include "match_consolidation.hpp"

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"DREAM-Stellar-match-consolidator", argc, argv};

    consolidation_arguments arguments;

    // Parser
    parser.info.author = "Evelin Aasna";
    parser.info.version = "1.0.0";

    parser.add_option(arguments.input_file,
                    'i',
                    "in",
                    "The file for GFF input.",
                    seqan3::option_spec::required,
                    seqan3::input_file_validator{{"gff"}});
    parser.add_option(arguments.output_file,
                    'o',
                    "out",
                    "The file for GFF output.",
                    seqan3::option_spec::required,
                    seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"gff"}});
    parser.add_option(arguments.overlap_length,
                    '\0',
                    "overlap",
                    "Segment overlap.",
                    seqan3::option_spec::required,
                    seqan3::arithmetic_range_validator{0, 500});
    parser.add_flag(arguments.verbose,
                    'v',
                    "verbose",
                    "Give more detailed information here.");

    try
    {
         parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "Parsing error. " << ext.what() << "\n";
        return -1;
    }

    consolidate_matches(arguments);

    return 0;
}
