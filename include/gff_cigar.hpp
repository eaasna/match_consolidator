#pragma once

#include "shared.hpp"

struct gff_cigar
{
    private:
        struct gff_cigar_operation
        {
            uint32_t pos;
            char op;

            gff_cigar_operation(uint32_t const p, char const o) : pos(p), op(o) {};

            void decrease_pos(size_t const decrease)
            {
                pos -= decrease;
            }

            std::string to_string() const
            {
                return std::to_string(pos) + op;
            }
        };

        std::vector<gff_cigar_operation> vec;
        std::string prefix = "cigar=";
        std::vector<char> valid_chars{'M', 'I', 'D'};
        uint64_t abs_pos;
        uint16_t len;

        void erase_el(bool const delete_front, size_t const el_to_delete)
        {
            if (delete_front)
            {
                for (size_t i = 0; i < el_to_delete; i++)
                    vec.erase(vec.begin());
            }
            else
            {
                for (size_t i = 0; i < el_to_delete; i++)
                    vec.pop_back();
            }
        }

    public:
        gff_cigar() = default;
        gff_cigar(gff_cigar const &) = default;
        gff_cigar & operator=(gff_cigar const &) = default;
        gff_cigar(gff_cigar &&) = default;
        gff_cigar & operator=(gff_cigar &&) = default;
        ~gff_cigar() = default;

        gff_cigar(uint64_t const dbegin, std::string const attribute_str, uint64_t const dend)
        {
            std::string cigar_str = attribute_str.substr(prefix.size());
            size_t first_ind = 0;
            for (size_t char_ind = 0; char_ind < cigar_str.size(); char_ind++)
            {
                if (std::count(valid_chars.begin(), valid_chars.end(), cigar_str[char_ind]) > 0)
                {
                    vec.emplace_back(static_cast<uint32_t>(std::stoi(cigar_str.substr(first_ind, char_ind))), cigar_str[char_ind]);
                    first_ind = char_ind + 1;
                }
            }

            abs_pos = dbegin;
            len = dend - dbegin + 1; // closed interval [dbegin; dend]
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
                size_t el_to_delete{0};
                if (vec[i].pos > shift)
                {
                    erase_el(delete_front, el_to_delete);
                    vec[i - el_to_delete] = gff_cigar_operation(vec[i].pos - static_cast<uint8_t>(shift), vec[i].op);
                    return;
                }
                else if (vec[i].pos == shift)
                {
                    el_to_delete++;
                    erase_el(delete_front, el_to_delete);
                    return;
                }
                else
                {
                    el_to_delete++;
                    shift -= vec[i].pos;
                }
            }

            return;
        }

        void join_cigars(gff_cigar const & other)
        {
            /*
            segments |------------|
                              |------------|
            this            xxxxxx
            other              xxx{xxx}

            join {xxx} to this
            */
            size_t join_from = abs_pos + len - other.abs_pos;
            size_t cum_sum{0}, ind{0};
            for (; ind < other.vec.size() && cum_sum < join_from; ind++)
            {
                cum_sum += static_cast<uint8_t>(other.vec[ind].pos);
            }

            if (vec.back().op == other.vec[ind-1].op)
                vec.back() = gff_cigar_operation(vec.back().pos + other.vec[ind-1].pos - join_from, vec.back().op);
            else
                vec.emplace_back(other.vec[ind-1].pos - join_from, other.vec[ind-1].op);

            vec.insert(vec.end(), other.vec.begin() + ind, other.vec.end());
        }

        std::string to_string()
        {
            std::string cigar_str = prefix;
            for (auto & c : vec)
            {
                cigar_str +=c.to_string();
            }

            return cigar_str;
        }
};
