#pragma once

#include <cmath>

#include "shared.hpp"

struct gff_percid
{
    //!TODO: if percid of joined alignment > allowed error rate -> can not join

    private:
        float perc;

    public:
        gff_percid() = default;
        gff_percid(gff_percid const &) = default;
        gff_percid & operator=(gff_percid const &) = default;
        gff_percid(gff_percid &&) = default;
        gff_percid & operator=(gff_percid &&) = default;
        ~gff_percid() = default;

        gff_percid(std::string const field5)
        {
            perc = std::stof(field5);
        }

        float get() const
        {
            return perc;
        }

        void update(uint64_t const dbegin, gff_mutations const & mut, uint64_t const dend)
        {
            // stellar coordinates are one-based closed intervals [dbegin, dend]
            // https://github.com/seqan/seqan/blob/f5f658343c366c9c3d44ba358ffc9317e78a09ed/apps/stellar/stellar_output.h#L192
            uint16_t ali_len = dend - dbegin + 1;
            perc = (1 - ( float ) mut.count() / ali_len) * 100;
        }

        std::string to_string() const
        {
            std::stringstream stream;
            // Stellar outputs floating point percentages with 4 decimal precision
            stream << std::fixed << std::setprecision(4) << perc;
            std::string str_percid = stream.str();
            // Remove trailing 0s
            str_percid.erase ( str_percid.find_last_not_of('0') + 1, std::string::npos );
            str_percid.erase ( str_percid.find_last_not_of('.') + 1, std::string::npos );

            return str_percid;
        }
};
