#pragma once

#include "shared.hpp"

struct gff_mutations
{
    gff_mutations() = default;
    gff_mutations(gff_mutations const &) = default;
    gff_mutations & operator=(gff_mutations const &) = default;
    gff_mutations(gff_mutations &&) = default;
    gff_mutations & operator=(gff_mutations &&) = default;
    ~gff_mutations() = default;

    struct mutation
    {
        uint16_t pos;
        char allele;

        mutation(std::string mutation_str)
        {
            pos = std::stoi(mutation_str.substr(0, mutation_str.size() - 1));
            allele = mutation_str[mutation_str.size() - 1];
        }

        void increase_pos(uint16_t const increase)
        {
            pos += increase;
        }
    };

    uint64_t abs_pos;
    uint16_t max_pos{0};
    std::vector<mutation> mut_vec;
    std::string prefix = "mutations=";

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
        auto mut_it = std::find_if(other.mut_vec.begin(),
                                other.mut_vec.end(),
                                [&](const auto & mut) { return mut.pos + other.abs_pos > abs_pos + max_pos; });

        for (;mut_it != other.mut_vec.end(); mut_it++)
        {
            mutation adjusted_mut = *mut_it;
            adjusted_mut.increase_pos(other.abs_pos - abs_pos);
            mut_vec.push_back(adjusted_mut);
        }
    }

    std::string to_string()
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
