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

#include "ibm.hpp"

#include <fstream>
#include <algorithm>
#include <limits>

#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/optimization/simplify.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/functions/ibm_helper.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/variable.hpp>
#include <cli/commands/ibm.hpp>
// #include <cli/commands/tvc.hpp>
#include <reversible/functions/ibm_helper.hpp>



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

ibm_command::ibm_command( const environment::ptr& env )
    : cirkit_command( env, "Translate Clifford+T circuits to IBM Q\nArchitecture: qx2 (default) or qx4" )
{
    opts.add_options()
    ( "all_perm,a",  "Try all permutations" )
    ( "rm_dup,r",  "Remove duplicate gates" )
    ( "ibm_qx4,4", "The IBM Qx4 is the target")
    ( "verbose,v",  "verbose" )
    ( "swap,s", "use swap based instead of template transformations")
    ( "matrix,p", "just print to matrix")
    ;
  add_new_option();
}


command::rules_t ibm_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

bool ibm_command::execute()
{
    std::vector<std::vector<unsigned>> qx2 ={{0,0,0,10,10},{4,0,0,10,10},{4,4,0,4,4},{10,10,0,0,0},{10,10,0,4,0}};
    std::vector<std::vector<unsigned>> qx4 ={{0,4,4,7,7},{0,0,4,7,7},{0,0,0,4,4},{3,3,0,0,0},{3,3,0,4,0}};

    auto& circuits = env->store<circuit>();
    circuit circ_working = circuits.current();
    circuit circ_IBM;
    // std::cout << " igates: " << circ_working.num_gates();
    unsigned start = circ_working.lines()+1;
    for(unsigned i = start ; i <= 5u; i++)
    {
        add_line_to_circuit( circ_working, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }
    
    if( is_set( "matrix") )
    {
        int circMatrix[5][5];
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                circMatrix[i][j] = 0;
        int maior = 0;
        for ( const auto& gate : circ_working ){
            if ( is_toffoli( gate ) && gate.controls().size() == 1 ){
                int linha = gate.controls().front().line();
                int coluna = gate.targets().front(); 
                circMatrix[linha][coluna]++;
                if (circMatrix[linha][coluna] > maior){
                    maior = circMatrix[linha][coluna];
                }
            }
        }

        int matrixVector[5];
        int maiorVector = 0;
        // for (int i = 0; i < 5; ++i)
        //     for (int j = 0; j < 5; ++j)
        //         std::cout << circMatrix[i][j] << " ";
        // std::cout << std::endl;

        for (int i = 0; i < 5; ++i)
        {
            matrixVector[i] = 0;
            for (int j = 0; j < 5; ++j)
            {
                matrixVector[i] = matrixVector[i] + circMatrix[i][j] + circMatrix[j][i];
            }
            if (matrixVector[i] > maiorVector)
                maiorVector = matrixVector[i];
        }
        // for (int i = 0; i < 5; ++i)
        //     for (int j = 0; j < 5; ++j)
        //         std::cout << float(float(circMatrix[i][j])/float(maior)) << " ";
        // std::cout << std::endl;
        for (int i = 0; i < 5; ++i)
            std::cout << matrixVector[i] << " ";
            // std::cout << float(float(matrixVector[i])/float(maiorVector)) << " ";
        // std::cout << std::endl;
    }

    if( !is_set( "all_perm" ) )
    {
        if ( is_set( "ibm_qx4" ) )
        {
            circ_working = transform_tof_clif(circ_working, qx4, 3);
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx4, !is_set( "swap" ) );
        }
        else
        {
            circ_working = transform_tof_clif(circ_working, qx2, 3);
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx2, !is_set( "swap" ) );
        }
        
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        if ( is_set( "rm_dup" ) )
        {
            circ_IBM = remove_dup_gates( circ_IBM);
        }
        circuits.current() = circ_IBM;
    }
    else // go through all permutations
    {
        int perm[5] = {0, 1, 2, 3, 4}, inv_perm[5], best_perm[5] = {0, 1, 2, 3, 4};
        unsigned best_cost = UINT_MAX;
        auto xx = circ_working.num_gates();
        circuit circ_best, permuted;
        std::vector<std::vector<int>> bp;
        do
        {
            clear_circuit(permuted);
            copy_circuit(circ_working, permuted);
            permute_lines( permuted, perm );
            if ( is_set( "ibm_qx4" ) )
            {
                permuted = transform_tof_clif(permuted, qx4, 3);
                circ_IBM = transform_to_IBMQ( permuted, map_method_qx4, !is_set( "swap" ) );
            }
            else
            {
                permuted = transform_tof_clif(permuted, qx2, 3);
                circ_IBM = transform_to_IBMQ( permuted, map_method_qx2, !is_set( "swap" ) );
            }
            if ( is_set( "new" ) )
            {
                circuits.extend();
            }
            
            if ( is_set( "rm_dup" ) )
            {
                circ_IBM = remove_dup_gates( circ_IBM );
            }
            circuits.current() = circ_IBM;
            if( is_set( "verbose" ) )
            {
                for( int i = 0; i < 5; i++ )
                {
                    std::cout << perm[i] << " ";
                }
                std::cout << "gates = " << circ_IBM.num_gates() << std::endl;
            }
            
            if( best_cost > circ_IBM.num_gates() )
            {
                bp.clear();
                // std::cout << "new best_cost = " << circ_IBM.num_gates() << "\n";
                best_cost = circ_IBM.num_gates();
                circ_best = circ_IBM;
                std::vector<int> v(5,0);
                for( int i = 0; i < 5; i++ )
                {
                    best_perm[i] = perm[i];
                    v[i] = perm[i];
                }
                bp.push_back(v);
            }
            else if( best_cost == circ_IBM.num_gates() )
            {
                std::vector<int> v(5,0);
                for( int i = 0; i < 5; i++ )
                {
                    v[i] = perm[i];
                }
                bp.push_back(v);
            }

            // // undo the permutation
            // for( int i = 0; i < 5; i++ )
            // {
            //     inv_perm[perm[i]] = i;
            // }
            // permute_lines( circ_working , inv_perm );
        } while ( std::next_permutation(perm,perm+5) );
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        circuits.current() = circ_best;
        std::cout << "best permutation = ";
        for( int i = 0; i < 5; i++ )
        {
            std::cout << best_perm[i] << " ";
        }
        std::cout << "gates = " << best_cost << std::endl;
        std::cout << "gates = " << best_cost - xx << std::endl;

        
        if( is_set( "matrix") )
        {
            int circMatrix[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    circMatrix[i][j] = 0;
            int maior = 0;
            for ( const auto& gate : circ_working ){
                if ( is_toffoli( gate ) && gate.controls().size() == 1 ){
                    int linha = gate.controls().front().line();
                    int coluna = gate.targets().front(); 
                    circMatrix[linha][coluna]++;
                    if (circMatrix[linha][coluna] > maior){
                        maior = circMatrix[linha][coluna];
                    }
                }
            }

            int matrixVector[5];
            int maiorVector = 0;
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    std::cout << circMatrix[i][j] << " ";
            for (int i = 0; i < 5; ++i){
                matrixVector[i] = 0;
                for (int j = 0; j < 5; ++j){
                    matrixVector[i] = matrixVector[i] + circMatrix[i][j] + circMatrix[j][i];
                }
                if (matrixVector[i] > maiorVector)
                    maiorVector = matrixVector[i];
            }
            for( int i = 0; i < 5; i++ )
                std::cout << best_perm[i] << " ";
            std::cout << std::endl;
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    std::cout << float(float(circMatrix[i][j])/float(maior)) << " ";
            std::cout << std::endl;
            for (int i = 0; i < 5; ++i)
                std::cout << float(float(matrixVector[i])/float(maiorVector)) << " ";
            std::cout << std::endl;
        }
    }
    // std::cout << "gates = " << circ_IBM.num_gates() << std::endl;
        
    return true;
}

command::log_opt_t ibm_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}



}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
