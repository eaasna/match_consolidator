#pragma once

#include "shared.hpp"

#include <seqan3/alphabet/cigar/cigar.hpp>

struct gff_cigar
{
    private:
        std::vector<seqan3::cigar> vec;
        std::string prefix = "cigar=";
        std::vector<char> valid_chars{'M', 'I', 'D'};

    public:
        gff_cigar() = default;
        gff_cigar(gff_cigar const &) = default;
        gff_cigar & operator=(gff_cigar const &) = default;
        gff_cigar(gff_cigar &&) = default;
        gff_cigar & operator=(gff_cigar &&) = default;
        ~gff_cigar() = default;

        gff_cigar(std::string const attribute_str)
        {
            using seqan3::get;
            using namespace seqan3::literals;

            std::string cigar_str = attribute_str.substr(prefix.size());
            size_t first_ind = 0;
            for (size_t char_ind = 0; char_ind < cigar_str.size(); char_ind++)
            {
                if (std::count(valid_chars.begin(), valid_chars.end(), cigar_str[char_ind]) > 0)
                {
                    seqan3::cigar::operation letter{};
                    letter.assign_char(cigar_str[char_ind]);
                    seqan3::cigar cig{stoi(cigar_str.substr(first_ind, char_ind)), letter};

                    vec.push_back(cig);
                    first_ind = char_ind + 1;
                }
            }
        }

        std::string to_string()
        {
            std::string cigar_str = prefix;
            for (auto & c : vec)
            {
                cigar_str += std::to_string(static_cast<uint32_t>(get<0>(c)));
                cigar_str += get<1>(c).to_char();
            }

            return cigar_str;
        }
};
