/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
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

#include "revtest.hpp"

#include <fstream>

#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <reversible/cli/stores.hpp>
#include <reversible/optimization/simplify.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/variable.hpp>

bool valid_CNOT[5][5] = {{0,1,1,0,0}, {0,0,1,0,0}, {0,0,0,0,0}, {0,0,1,0,1}, {0,0,1,0,0}};
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

revtest_command::revtest_command( const environment::ptr& env )
  : cirkit_command( env, "Reversible circuit simplification" )
{
  opts.add_options()
    ( "methods",   value_with_default( &methods ), "optimization methods:\nm: try to merge gates with same target\nn: cancel NOT gates\na: merge adjacent gates\ne: resynthesize same-target gates with exorcism\ns: propagate SWAP gates (may change output order)" )
    ( "noreverse",                                 "do not optimize in reverse direction" )
    ;
  be_verbose();
  add_new_option();
}


command::rules_t revtest_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

    
bool revtest_command::execute()
{
    auto& circ = env->store<circuit>();

//    const auto vars = create_name_list( "v%d", circ.lines() );
    const auto vars = create_name_list( "v%d", 5u );
    
    // iterate through the gates
    for ( const auto& gate : *circ )
    {
        unsigned target = gate.targets().front();
        if ( is_toffoli( gate ) )
        {
            if( gate.controls().empty() ) // a NOT gate
            {
                std::cout << "x q[" << target << "];" << std::endl;
            }
            else // CNOT gate
            {
                unsigned control = gate.controls().front().line();
                
                if( valid_CNOT[control][target] )
                {
                    std::cout << "cx q[" << control << "],q[" << target << "];" << std::endl;
                }
                else if( valid_CNOT[target][control] ) // invert CNOT
                {
                    std::cout << "h q[" << control << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                    std::cout << "cx q[" << target << "],q[" << control << "];" << std::endl;
                    std::cout << "h q[" << control << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                }
                else // swap target with 2
                {
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    std::cout << "h q[" << 2 << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    std::cout << "h q[" << 2 << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    
                    std::cout << "cx q[" << control << "],q[" << 2 << "];" << std::endl;
                    
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    std::cout << "h q[" << 2 << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    std::cout << "h q[" << 2 << "];"  << std::endl;
                    std::cout << "h q[" << target << "];"  << std::endl;
                    std::cout << "cx q[" << target << "],q[" << 2 << "];" << std::endl;
                    
                }
            }
        }
        else if ( is_pauli( gate ) )
        {
            const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            
            switch ( tag.axis )
            {
                case pauli_axis::X:
                    assert( tag.root == 1u );
                    std::cout << "x";
                    break;
                    
                case pauli_axis::Y:
                    assert( tag.root == 1u );
                     std::cout << "y";
                    break;
                    
                case pauli_axis::Z:
                    switch ( tag.root )
                {
                    case 1u:
                        std::cout << "z";
                        break;
                    case 2u:
                        std::cout << "s";
                        break;
                    case 4u:
                        std::cout << "t";
                        break;
                    default:
                        assert( false );
                }
                    break;
            }
            
            if ( tag.adjoint )
            {
                std::cout << "dg";
            }
            
             std::cout << " q[" << target << "];" << std::endl;
        }
        else if ( is_hadamard( gate ) )
        {
             std::cout << "h q[" << target << "];" << std::endl;
        }
        else
        {
            assert( false );
        }
    }
    
    
    const auto& circuits = env->store<circuit>();
    circuit circ_IBM = transform_to_IBM_Q5( circuits.current() );
    
    write_qc( circ_IBM, "../../test.qc", false );
    return true;
}

command::log_opt_t revtest_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

// transform a Clifford+T circuit to be IBM compliant
circuit transform_to_IBM_Q5( const circuit& circ )
{
    circuit circ_IBM;
    unsigned target, control;
    std::vector<unsigned int> new_controls;
    copy_metadata( circ, circ_IBM );
    
    // iterate through the gates
    for ( const auto& gate : circ )
    {
        
        target = gate.targets().front();
        new_controls.clear();
        new_controls.push_back( target );
        if( !gate.controls().empty() )
        {
            control = gate.controls().front().line();
        }
        
        if ( is_toffoli( gate ) )
        {
            if( gate.controls().empty() ) // a NOT gate
            {
                std::cout << "NOT ";
                append_toffoli( circ_IBM, gate.controls(), target );
            }
            else // CNOT gate
            {
                
                if( valid_CNOT[control][target] )
                {
                    append_toffoli( circ_IBM, gate.controls(), target );
                }
                else if( valid_CNOT[target][control] ) // invert CNOT
                {
                    append_hadamard( circ_IBM, control );
                    append_hadamard( circ_IBM, target );
                    append_toffoli( circ_IBM, new_controls, control );
                    append_hadamard( circ_IBM, control );
                    append_hadamard( circ_IBM, target );
                                    }
                else // swap target with 2
                {
                    append_toffoli( circ_IBM, new_controls, 2u );
                    append_hadamard( circ_IBM, 2u );
                    append_hadamard( circ_IBM, target );
                    append_toffoli( circ_IBM, new_controls, 2u );
                    append_hadamard( circ_IBM, 2u );
                    append_hadamard( circ_IBM, target );
                    append_toffoli( circ_IBM, new_controls, 2u );
                    
                    append_toffoli( circ_IBM, gate.controls(), 2u );
                    
                    append_toffoli( circ_IBM, new_controls, 2u );
                    append_hadamard( circ_IBM, 2u );
                    append_hadamard( circ_IBM, target );
                    append_toffoli( circ_IBM, new_controls, 2u );
                    append_hadamard( circ_IBM, 2u );
                    append_hadamard( circ_IBM, target );
                    append_toffoli( circ_IBM, new_controls, 2u );
                    
                }
            }
        }
        else if ( is_pauli( gate ) )
        {
            
            const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            
            append_pauli( circ_IBM, target, tag.axis, tag.root, tag.adjoint );
 /*           
            switch ( tag.axis )
            {
                case pauli_axis::X:
                    assert( tag.root == 1u );
                    std::cout << "x";
                    break;
                    
                case pauli_axis::Y:
                    assert( tag.root == 1u );
                    std::cout << "y";
                    break;
                    
                case pauli_axis::Z:
                    switch ( tag.root )
                {
                    case 1u:
                        std::cout << "z";
                        break;
                    case 2u:
                        std::cout << "s";
                        break;
                    case 4u:
                        std::cout << "t";
                        break;
                    default:
                        assert( false );
                }
                    break;
            }

            if ( tag.adjoint )
            {
                std::cout << "dg";
            }
            
            std::cout << " q[" << target << "];" << std::endl;
 */
        }
        else if ( is_hadamard( gate ) )
        {
            append_hadamard( circ_IBM, target );
        }
        else
        {
            assert( false );
        }
    }

    
    return circ_IBM;
}
    
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End: