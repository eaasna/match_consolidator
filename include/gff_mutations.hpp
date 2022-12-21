#pragma once

#include "shared.hpp"

struct gff_mutations
{
    private:
        struct mutation
        {
            uint16_t pos;
            char allele;

            mutation(std::string mutation_str)
            {
                pos = std::stoi(mutation_str.substr(0, mutation_str.size() - 1));
                allele = mutation_str[mutation_str.size() - 1];
            }

            template <typename int_t>
            void adjust_pos(int_t const adjustment)
            {
                pos += adjustment;
            }
        };

        std::vector<mutation> mut_vec;
        uint64_t abs_pos;
        uint16_t max_pos{0};
        std::string prefix = "mutations=";

    public:
        gff_mutations() = default;
        gff_mutations(gff_mutations const &) = default;
        gff_mutations & operator=(gff_mutations const &) = default;
        gff_mutations(gff_mutations &&) = default;
        gff_mutations & operator=(gff_mutations &&) = default;
        ~gff_mutations() = default;

        gff_mutations(std::string const mutations_str, uint64_t const dbegin)
        {
            std::istringstream iss(mutations_str.substr(prefix.size()));
            std::string field;
            while(std::getline(iss, field, ','))
            {
                mutation mut(field);
                mut_vec.emplace_back(mut);
                if (mut.pos > max_pos)
                    max_pos = mut.pos;
            }

            abs_pos = dbegin;
        }

        void join_mutations(gff_mutations const & other)
        {
            /*
            segments |------------|
                              |------------|
            this            xxxxxx
            other              xxx{xxx}

            join {xxx} to this
            */
            auto mut_it = std::find_if(other.mut_vec.begin(),
                                    other.mut_vec.end(),
                                    [&](const auto & mut) { return mut.pos + other.abs_pos > abs_pos + max_pos; });

            for (;mut_it != other.mut_vec.end(); mut_it++)
            {
                mutation adjusted_mut = *mut_it;
                adjusted_mut.adjust_pos<uint16_t>(other.abs_pos - abs_pos);
                mut_vec.push_back(adjusted_mut);
            }
        }

        size_t count() const
        {
            return mut_vec.size();
        }

        std::pair<bool, uint16_t> shorten_match_greedily(uint64_t const dend)
        {
            uint16_t begin_shift = mut_vec[0].pos;
            uint16_t end_shift = dend + 1 - max_pos;

            if (begin_shift < end_shift)
            {
                abs_pos += begin_shift + 1;
                mut_vec.erase(mut_vec.begin());

                for (auto & mut : mut_vec)
                    mut.adjust_pos<int16_t>(begin_shift * (-1));

                return std::make_pair(true, begin_shift);
            }
            else
            {
                mut_vec.erase(mut_vec.end() - 1);
                max_pos = mut_vec[-1].pos;

                return std::make_pair(false, end_shift);
            }
        }

        std::string to_string() const
        {
            std::string mut_str = prefix;

            for (auto mut : mut_vec)
            {
                mut_str += std::to_string(mut.pos);
                mut_str += mut.allele;
                mut_str += ",";
            }

            return mut_str.substr(0, mut_str.size() - 1);
        }
};
