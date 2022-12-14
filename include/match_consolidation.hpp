#pragma once

#include "shared.hpp"

#include <seqan3/io/sequence_file/all.hpp>

/*! \brief Function, consolidating matches after distributed Stellar alignment.
 *
 *  TODO:
 *  Pick matches based on search scheme (single-best, all-best etc)
 */
void consolidate_matches(consolidation_arguments const & arguments);
