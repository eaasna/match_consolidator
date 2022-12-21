#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "gff_mutations.hpp"
#include "gff_percid.hpp"

TEST(too_many_errors, ok)
{
    gff_percid percid("99.0", 0.05);

    EXPECT_FALSE(percid.too_many_errors());
}

TEST(too_many_errors, edge_case)
{
    gff_percid percid("99.0", 0.01);

    EXPECT_TRUE(percid.too_many_errors());
}

TEST(too_many_errors, too_many)
{
    gff_percid percid("99.0", 0.0);

    EXPECT_TRUE(percid.too_many_errors());
}
