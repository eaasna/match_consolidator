#pragma once

#include "shared.hpp"

#include <seqan3/alphabet/cigar/cigar.hpp>

struct gff_cigar
{
    private:
        std::vector<seqan3::cigar> vec;
        std::string prefix = "cigar=";
        std::vector<char> valid_chars{'M', 'I', 'D'};

        void erase_el(bool const delete_front)
        {
            if (delete_front)
                vec.erase(vec.begin());
            else
                vec.pop_back();
        }

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
                    seqan3::cigar cig{static_cast<uint32_t>(std::stoi(cigar_str.substr(first_ind, char_ind))), letter};
                    vec.push_back(cig);
                    first_ind = char_ind + 1;
                }
            }
        }

        void adjust_pos(bool const delete_front, size_t shift)
        {
            int8_t incr = 1;
            int8_t i = 0;
            if (!delete_front)
            {
                incr = -1;
                i = vec.size() - 1;
            }

            for (;i < static_cast<uint16_t>(vec.size()) && i > -1; i += incr)
            {
                if (get<0>(vec[i]) > shift)
                {
                    seqan3::cigar new_cig{get<0>(vec[i]) - static_cast<uint8_t>(shift), get<1>(vec[i])};
                    vec[i] = new_cig;
                    return;
                }
                else if (get<0>(vec[i]) == shift)
                {
                    erase_el(delete_front);
                    return;
                }
                else
                {
                    erase_el(delete_front);
                    shift -= get<0>(vec[i]);
                }
            }

            return;
        }

        void join_cigars(gff_cigar const & other)
        {

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
