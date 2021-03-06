/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "rms.hpp"

#include <alice/rules.hpp>
#include <core/utils/program_options.hpp>

#include <cli/reversible_stores.hpp>
#include <reversible/synthesis/reed_muller_synthesis.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

rms_command::rms_command( const environment::ptr& env )
  : cirkit_command( env, "Reed-Muller based synthesis" )
{
  opts.add_options()
    ( "bidirectional,b", value_with_default( &bidirectional ), "bidirectional synthesis" )
    ;
  add_new_option();
  be_verbose();
}

command::rules_t rms_command::validity_rules() const
{
  return {has_store_element<binary_truth_table>( env )};
}

bool rms_command::execute()
{
  const auto& specs = env->store<binary_truth_table>();
  auto& circuits = env->store<circuit>();

  extend_if_new( circuits );

  auto settings = make_settings();
  settings->set( "bidirectional", bidirectional );
  reed_muller_synthesis( circuits.current(), specs.current(), settings, statistics );

  print_runtime();

  return true;
}

command::log_opt_t rms_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
