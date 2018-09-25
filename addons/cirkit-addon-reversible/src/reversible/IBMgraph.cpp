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

#include "IBMgraph.hpp"
#include <reversible/io/print_circuit.hpp>
#include <reversible/functions/move_qubit.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <boost/lexical_cast.hpp>
#include <reversible/functions/ibm_helper.hpp>


namespace cirkit
{
    using matrix = std::vector<std::vector<int>>;
    // Here the memory error stopped
    std::vector<TransPath> path_list;
    bool **graph_adjacency = NULL;  // graph structure
    int **trans_cost = NULL;        // the cost of each cnot
    TransPath **trans_path = NULL;  // the transformation path or a given cnot
    
    bool read_graph( const std::string& filename )
    {
        int v,w;
        std::fstream graphfile;
        graphfile.open ( filename, std::fstream::in );
        if ( !graphfile.is_open() )
        {
            return false;
        }
        graphfile >> graph_size;
        graph_adjacency = new bool*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            graph_adjacency[i] = new bool[graph_size];
            for(int j = 0; j < graph_size; j++)
            {
                graph_adjacency[i][j] = false;
            }
        }
        while ( !graphfile.eof() )
        {
            graphfile >> v >> w;
            graph_adjacency[v][w] = true;
        }
        graphfile.close();
        return true;
    }

    // read the matrix and the transformations from a file
    bool read_from_file ( const std::string& filename )
    { 
        std::vector<std::string> type_name = { "cab", "cba", "tab", "tba", "cabi", "cbai", "tabi", "tbai", "nop", "flip", "cnot3"};
        std::string tmp;
        TransPath tp;
        int a,b,c,v,w;
        std::ifstream graphfile ( filename );
        if ( !graphfile.is_open() )
            return false;
        
        graphfile >> graph_size; // get the number of qubits
        graph_adjacency = new bool*[graph_size];
        allocate_data_stuctures();

        // read the graph
        for(int i = 0; i < graph_size ; i++ ){
            graph_adjacency[i] = new bool[graph_size];
            for(int j = 0; j < graph_size; j++){
                graphfile >> tmp;
                if(tmp == "X")
                    graph_adjacency[i][j] = true;
                else
                    graph_adjacency[i][j] = false;    
            }
        }
        // read the matrix with the transformations costs
        for( v = 0; v < graph_size; v++){ 
            for( w = 0; w < graph_size; w++)
                graphfile >> trans_cost[v][w];
        }
        
        v = 0; // for the variable trans_path
        w = 1; // starting with cnot(0,1)
        while ( !graphfile.eof() ){
            graphfile >> tmp;
            if(tmp == "cost"){ // it means that the transformation has already been read
                trans_path[v][w] = tp;  // update trans_path,
                tp.clear();             // and clear for the next one
                ++w;
                if(v == w) // cnot(x,x) doesn't exist
                    ++w;
                if(w == graph_size){    // finished all w possibilities
                    ++v;                // and start the next one
                    w = 0;
                }
            }
            auto it = std::find(type_name.begin(), type_name.end(), tmp);   
            if(it != type_name.end()){                                  // check if it is a movement
                unsigned pos = std::distance(type_name.begin(), it);    // get the type of movement
                if(tmp == "cnot3"){     // cnot3 needs three parameters
                    graphfile >> a >> b >> c;
                    tp.add( MoveQubit( pos, a, b, c ));
                }
                else{                   // the other movements only two
                    graphfile >> a >> b;
                    tp.add( MoveQubit( pos, a, b ));
                }
            }
        }
        graphfile.close();
        return true;
    }

    // write the graph, the matrix and the transformations to a file
    bool write_to_file ( const std::string& filename )
    {
        std::ofstream graphfile ( filename );
        if ( !graphfile.is_open() )
            return false;

        graphfile << graph_size << std::endl;

        // writing the graph
        for( int i = 0; i < graph_size; i++ ){ 
            for( int j = 0; j < graph_size; j++ ){
                graphfile<< (graph_adjacency[i][j] ? "X " : "- ");
            }
            graphfile << std::endl;
        }
        // writing the matrix
        for( int v = 0; v < graph_size; v++){   
            for( int w = 0; w < graph_size; w++)
                graphfile << trans_cost[v][w] << " ";
            graphfile << std::endl;
        }
        // writing the transformations
        for( int v = 0; v < graph_size; v++){
            for( int w = 0; w < graph_size; w++){
                if( v != w ){
                    graphfile << "cnot(" << v << "," << w << ") => ";
                    trans_path[v][w].print( graphfile );
                }
            }
        }
        graphfile.close();
        return true;
    }

   void print_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            for( int j = 0; j < graph_size; j++ ){
                std::cout << (graph_adjacency[i][j] ? "X " : "- ");
            }
            std::cout << std::endl;
        }
    }

    // print the matrix and the transformations
    void print_matrix( )
    {
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
                std::cout << trans_cost[v][w] << " ";
            std::cout << std::endl;
        }
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                if( v != w )
                {
                    std::cout << "cnot(" << v << "," << w << ") => ";
                    trans_path[v][w].print( );
                }
            }
        }
    }

    void delete_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            delete [] graph_adjacency[i];
            delete [] trans_path[i];
            delete [] trans_cost[i];
            
        }
        delete [] graph_adjacency;
        delete [] trans_path;
        delete [] trans_cost;
        graph_size = 0;
    }
    
    // When a path is found add it to the global vector of paths.
    void find_all_paths( int v, int w, TransPath &tp, bool *visited )
    {
        bool done = false;
        for( int i = 0; i < graph_size; i++ )
        {
            if( graph_adjacency[v][w] )
            {
                tp.add( MoveQubit( nop, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( graph_adjacency[w][v] )
            {
                tp.add( MoveQubit( flip, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( done ) return; // no need to go further
            if( graph_adjacency[v][i] && !visited[i] )
            {
                tp.add( MoveQubit( cab, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][v] && !visited[i] )
            {
                tp.add( MoveQubit( cba, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[w][i] && !visited[i] )
            {
                tp.add( MoveQubit( tab, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );           
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][w] && !visited[i] )
            {
                tp.add( MoveQubit( tba, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
        }
    }

    /* Precondition: the path_list contains all paths from v to w
     Postcondition: the best path and its cost will be stored in:
        - trans_cost[v][w]
        - trans_path[v][w]
     */
    
    void set_best_path(int v, int w )
    {
        TransPath tp, best_tp;
        int best_cost;
        best_tp = path_list[0];
        best_cost = best_tp.costPlus();

        for ( auto &p : path_list )
        {
            p.movCnot3();
            if(p.costPlus() < best_cost )
            {
                best_cost = p.costPlus();
                best_tp = p;
            }
        }

        trans_cost[v][w] = best_cost;
        best_tp.addInverse();
        trans_path[v][w] = best_tp;
    }

    void allocate_data_stuctures(){
        trans_cost = new int*[graph_size];
        trans_path = new TransPath*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            trans_cost[i] = new int[graph_size];
            trans_path[i] = new TransPath[graph_size];
        }
    }
    
    void create_trans( bool verbose )
    {
        allocate_data_stuctures();
        path_list.clear();
        bool visited[graph_size];
        for( int i = 0; i < graph_size ; i++ )
        {
            visited[i] = false;
        }

        TransPath tp;
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                if( v == w )
                {
                    trans_cost[v][w] = 0;
                }
                else{
                    for( int i = 0; i < graph_size ; i++ ) visited[i] = false;
                    visited[v] = true;
                    visited[w] = true;
                    path_list.clear();
                    tp.clear();
                    find_all_paths( v,  w, tp, visited );
                    set_best_path( v,  w);
                }
            }
        }
        if ( verbose )
        {
            for( int v = 0; v < graph_size; v++)
            {
                for( int w = 0; w < graph_size; w++)
                {
                    std::cout << trans_cost[v][w] << " ";
                }
                std::cout << std::endl;
            }
            // std::cout << "== Optimization ==" << std::endl;
            // for( int v = 0; v < graph_size; v++)
            // {
            //     for( int w = 0; w < graph_size; w++)
            //     {
            //         std::cout << trans_cost[v][w]-trans_path[v][w].opt() << " ";
            //     }
            //     std::cout << std::endl;
            // }
            for( int v = 0; v < graph_size; v++)
            {
                for( int w = 0; w < graph_size; w++)
                {
                    if( v != w )
                    {
                        std::cout << "cnot(" << v << "," << w << ") => ";
                        trans_path[v][w].print();
                        // std::cout << "Can reduce: " << trans_path[v][w].opt() << std::endl;
                    }
                }
            }
        }
    }

    // expand the cnot gates that are not supported by the architecture
    // assume that the corresponding matricies have been set up correctly
    void expand_cnots( circuit& circ_out, const circuit& circ_in ){
        
        unsigned target, control, moreCnot3 = 0;
        std::vector<unsigned int> new_controls, control2, old_controls;
        
        copy_metadata( circ_in, circ_out );
        for ( const auto& gate : circ_in )
        {
            target = gate.targets().front();
            new_controls.clear();
            new_controls.push_back( target );
            if( !gate.controls().empty() )
            {
                control = gate.controls().front().line();
                old_controls.clear();
                old_controls.push_back( control );
            }
            
            if ( is_toffoli( gate ) )
            {
                if( gate.controls().empty() ) // a NOT gate
                {
                    circ_out.append_gate() = gate;
                }
                else // CNOT gate
                {
                    moreCnot3 = 0;
                    for ( auto &p : trans_path[control][target].tpath ) {
                        switch ( p.getType() )
                        {
                            case cab : //std::cout << "cab" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case cba : //std::cout << "cba" << std::endl;
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case tab : //std::cout << "tab" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case tba : //std::cout << "tba" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                break;
                            case cabi : //std::cout << "cabi" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case cbai : //std::cout << "cbai" << std::endl;
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                break;
                            case tabi : //std::cout << "tabi" << std::endl;
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case tbai : //std::cout << "tbai" << std::endl;
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case nop : //std::cout << "nop" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case flip : //std::cout << "flip" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case cnot3 :
                                if(moreCnot3 == 0)  // if it is the first cnot3
                                {                   // just append the four cnot gates
                                    append_cnot( circ_out, p.getA(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getC() );
                                    append_cnot( circ_out, p.getA(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getC() );
                                    ++moreCnot3;    // update the number of cnot3
                                }
                                else    // if it is more than "two arrows"
                                {       // we have to calculate the number of cnots
                                    unsigned c = pow(2,moreCnot3) + pow(2,++moreCnot3) - 2; // the calculation: 2^n + 2^(n+1) - 2 -> n=number of "arrows" - 1
                                    append_cnot( circ_out, p.getB(), p.getC() );                    // append the cnot   
                                    for (int i = 0, j = circ_out.num_gates()-(c+1); i < c; ++i, ++j)// and copy all the cnots placed before
                                        circ_out.append_gate() = circ_out[j];
                                    append_cnot( circ_out, p.getB(), p.getC() );                    // append again the cnot
                                    ++moreCnot3;
                                }
                                break;
                            default : std::cout << "ERROR expand_cnots" << std::endl;
                        }
                    }
                }
            }
            else
            {
               circ_out.append_gate() = gate;
            }
        }
    }

    




void pretty_print(const std::vector<int>& v) {
  static int count = 0;
  std::cout << "combination no " << (++count) << ": [ ";
  for (int i = 0; i < v.size(); ++i) { std::cout << v[i] << " "; }
  std::cout << "] " << std::endl;
}

    unsigned custo;
    matrix mapeamento;
    matrix matrix_custo;

    void extract_matrix(std::vector<int>& v)
    {
        for (int i = 0; i < v.size(); ++i)
        {
            for (int j = 0; j < v.size(); ++j)
            {
                matrix_custo[i][j] = trans_cost[v[i]][v[j]];
                // std::cout << " " << matrix_q[i][j];
            }
            // std::cout << std::endl;
        }
    }

    void all_combinations(int offset, int k, std::vector<int>& qubits, std::vector<int>& combination, matrix& matrix_circuito) 
    {
        // long long int total = 0;
        
        if (k == 0) 
        {
            // pretty_print(combination);
            unsigned int aux;
            do
            {
                // ++total;
                // if(total%500000000 == 0)
                //     std::cout << " | " << std::endl;
                // std::cout << contagem++ << ": "; 
                // for (int i = 0; i < 5; ++i)
                // {
                //     std::cout << " " << v[i];
                // }
                // std::cout << std::endl;
                extract_matrix(combination);
                aux = 0;
                for (int i = 0; i < combination.size(); ++i)
                {
                    for (int j = 0; j < combination.size(); ++j)
                    {
                        aux += matrix_custo[i][j]*matrix_circuito[i][j];
                    }
                }
                if(aux < custo)
                {
                    custo = aux;
                    mapeamento.clear();
                    mapeamento.push_back(combination);
                }
                else if( aux == custo)
                {
                    mapeamento.push_back(combination);
                }
                // if(aux == 22)
                //     getchar();
                // custo.push_back(aux);
                // mapeamento.push_back(v);
                // if( aux < custo )
                // {
                //     custo = aux;
                //     for (int i = 0; i < 4; ++i)
                //     {
                //         std::cout << " " << v[i];
                //     }
                // }

            } while ( std::next_permutation(combination.begin(),combination.begin()+combination.size()) );
            return;
        }
        for (int i = offset; i <= qubits.size() - k; ++i) 
        {
            combination.push_back(qubits[i]);
            all_combinations(i+1, k-1, qubits, combination, matrix_circuito);
            combination.pop_back();
        }
    }

    void the_mapping( circuit& circ_out, const circuit& circ_in )
    {
        std::cout << "THE MAPPING!" << std::endl;
        circuit aux;
        copy_circuit(circ_in, aux);
        // for( int v = 0; v < graph_size; v++)
        // {
        //     for( int w = 0; w < graph_size; w++)
        //         std::cout << trans_cost[v][w] << " ";
        //     std::cout << std::endl;
        // }
        unsigned target, control;
        std::vector<int> qubits;
        std::vector<int> combination;
        matrix matrix_circuito;
        custo = 100000;


        for (int i = 0; i < mapeamento.size(); ++i)
            mapeamento[i].clear();
        mapeamento.clear();
        for (int i = 0; i < matrix_custo.size(); ++i)
            matrix_custo[i].clear();
        matrix_custo.clear();


        for (int i = 0; i < circ_in.lines(); ++i)
            combination.push_back(0);
        for (int i = 0; i < circ_in.lines(); ++i)
        {
            matrix_circuito.push_back(combination);
            matrix_custo.push_back(combination);
        }
        
        combination.clear();

        for ( const auto& gate : circ_in )
        {
            if( !gate.controls().empty() ) // if is not a NOT gate
            {
                target = gate.targets().front();
                control = gate.controls().front().line();
                ++matrix_circuito[control][target];
            }
        }

        for (int i = 0; i < graph_size; ++i) 
            qubits.push_back(i);

        all_combinations(0, circ_in.lines(), qubits, combination, matrix_circuito);

        // for (int i = 0; i < custo.size(); ++i) { 
        //     std::cout << custo[i] << std::endl;
        //     for (int j = 0; j < 5; ++j) {
        //         std::cout << " " << mapeamento[i][j];
        //     } 
        //     std::cout << std::endl;

        // }
        
        std::cout << "\nMenor: " << custo << std::endl;
        for (int i = 0; i < mapeamento.size(); ++i) 
        { 
            for (int j = 0; j < circ_in.lines(); ++j) 
            {
                std::cout << " " << mapeamento[i][j];
            } 
            std::cout << std::endl;
        }

        for(unsigned i = circ_in.lines() ; i < graph_size; i++)
            add_line_to_circuit( aux, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
        
        permute_lines(aux, &mapeamento[0][0]);

        // copy_circuit(circ_out, aux);
        expand_cnots( circ_out, aux );

        // std::cout << "Menor: " << *min_element(std::begin(custo), std::end(custo)) << std::endl;
        // std::vector<int>::iterator it = std::find(custo.begin(), custo.end(), 25);
        // int index = std::distance(custo.begin(), it);
        
        // for (int j = 0; j < 5; ++j) {
        //     std::cout << " " << mapeamento[index][j];
        // }
        // std::cout << std::endl;
        
        // std::cout << "Total: " << total << std::endl;

        // std::cout << "Max: " << *max_element(std::begin(custo), std::end(custo)) << std::endl;
        //  it = std::find(custo.begin(), custo.end(), 535);
        //  index = std::distance(custo.begin(), it);
        
        // for (int j = 0; j < 5; ++j) {
        //     std::cout << " " << mapeamento[index][j];
        // }
        // std::cout << std::endl;


        // std::vector<unsigned> s; // cnot implementable
        // std::vector<unsigned> r; // reverse of s
        // std::vector<unsigned> t; // total
        // std::pair<int,int> ncontrol;
        // std::pair<int,int> ntarget;

        // for(int i = 0; i < graph_size ; i++ )
        // {
        //     unsigned ns = 0;
        //     unsigned nr = 0;
        //     for(int j = 0; j < graph_size; j++)
        //     {
        //         if(graph_adjacency[i][j])
        //             ++ns;
        //         if(graph_adjacency[j][i])
        //             ++nr;
        //     }
        //     s.push_back(ns);
        //     r.push_back(nr);
        //     t.push_back(ns+nr);
        //     std::cout << "Qubit " << i << " - " << ns << " X's in this line and " << nr << " X's in this column"<< std::endl;
        // }
        // for (int i = 0; i < graph_size; ++i)
        //     std::cout << "Qubit " << i << ": " << t[i] << std::endl;
    }

    
}
