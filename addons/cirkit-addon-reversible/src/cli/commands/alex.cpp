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

#include "alex.hpp"

#include <cmath>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <iostream>

#include <boost/program_options.hpp>

#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include <cli/reversible_stores.hpp>
#include <reversible/utils/matrix_utils.hpp>
#include <alice/rules.hpp>


namespace cirkit
{

using boost::program_options::value;

alex_command::alex_command( const environment::ptr& env )
	: cirkit_command( env, "Alex test" )
{
	opts.add_options()
	( "filename,f",    value( &filename ),    "name of the output file" )
	;

}

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}


using matrix = std::vector<std::vector<unsigned>>;
std::ofstream outputFile;

// Return the higher value element in the matrix
int get_max_element(matrix& m, unsigned& l, unsigned& c)
{
    unsigned h = 0;
    for (int i = 0; i < m.size(); ++i)
    {
        for (int j = 0; j < m[i].size(); ++j)
        {
            if(i != j && m[i][j] > h)
            {
                h = m[i][j];
                l = i;
                c = j;
            }
        }
    }
    return h;
}

// return the number of different gates of the circuit
int getNumberDifGates( matrix& c )
{
	unsigned qtd = 0;
	for (int i = 0; i < c.size(); ++i)
		for (int j = 0; j < c[i].size(); ++j)
			if( c[i][j] > 0 )
				++qtd;
	return qtd;
}

// Function to print the objective function
void printObjectiveFunction( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	outputFile << "/* Begin Objective Function */" << std::endl;
	outputFile << "min:\t";
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							outputFile << qx[k][m]*cnot[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << ";";
						else
							outputFile << qx[k][m]*cnot[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				outputFile << std::endl << "\t";
				++aux;
			}
			else if(line)
			{
				outputFile << std::endl;
			}
		}
	}
	outputFile << "/* End Objective Function */" << std::endl;
}

// Function to print the first restriction
void printFirstRestriction( matrix& cnot )
{
	unsigned aux = 0;
	outputFile << "/* Begin First Restriction */" << std::endl;
	for (int i = 0; i < cnot.size(); ++i)
	{
		for (int j = 0; j < cnot.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < cnot.size(); ++k)
			{
				for (int m = 0; m < cnot.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == cnot.size()-1 && m == cnot.size()-2 )
							outputFile << cnot[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << " = " << cnot[i][j] << ";";
						else
							outputFile << cnot[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line)
				outputFile << std::endl;
		}
	}
	outputFile << "/* End First Restriction */" << std::endl;
}

// Function to print the final restriction
void printEndRestriction( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	outputFile << "/* Begin Final Restriction */" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " = " << difGates << ";";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line)
			{
				outputFile << std::endl;
				++aux;
			}
		}
	}
	outputFile << "/* End Final Restriction */" << std::endl;
}

// Function to print the variables
void printIntegerVariables( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	outputFile << "/* Begin Integer Variables */" << std::endl;
	outputFile << "int\t";
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << ";";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << "  ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				outputFile << std::endl << "\t";
				++aux;
			}
			else if(line)
			{
				outputFile << std::endl;
			}
		}
	}
	outputFile << "/* End Integer Variables */" << std::endl;
}

// Create a matrix with 0's
void createMatrix( matrix& m, unsigned size )
{
	// std::cout << "Creating matrix..." << std::endl;
  	std::vector<unsigned> v;
	for (int i = 0; i < size; ++i)
		v.push_back(0);
	for (int i = 0; i < size; ++i)
		m.push_back(v);
}

// Create a matrix with the cnots 
void generateMatrixCnots( circuit& circ, matrix& m )
{
	// std::cout << "Generating matrix..." << std::endl;	
  	unsigned target, control;
	for ( const auto& gate : circ )
	{
		if( !gate.controls().empty() )
		{
		  target = gate.targets().front();
		  control = gate.controls().front().line();
		  ++m[control][target];
		}
	}
}

// Print the cnots in the circuit
void printMatrixCnots( matrix& m )
{
	// std::cout << "Printing matrix..." << std::endl;
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
			std::cout << "\t" << m[i][j];
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// Only function to write almost everything
// type = 0 -> control - control
// type = 1 -> target - target
// type = 2 -> control - target
// type = 3 -> target - control
// type = 4 -> inverse (control -> target; target -> control)
// type = 13 -> Ghost gates
void writeDep( unsigned l0, unsigned c0, unsigned l1, unsigned c1, unsigned size, unsigned type )
{
	switch ( type )
    {
    	case 0:
			outputFile << "/* Writing Control - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 1:
			outputFile << "/* Writing Target - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 2:
			outputFile << "/* Writing Control - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 3:
			outputFile << "/* Writing Target - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 4:
			outputFile << "/* Writing Inverse dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	default:
			outputFile << "/* ERRORRRRRRRRRRRRRRRRR */" << std::endl;
    }
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			if(i != j && type != 4)
			{
				if(type == 0 || type == 2)
				{
					outputFile << "G" << l0 << "_" << c0;
					outputFile << "c" << i << "_" << j;
					outputFile << " <= ";
				}
				else
				{
					outputFile << "G" << l0 << "_" << c0;
					outputFile << "c" << j << "_" << i;
					outputFile << " <= ";
				}
				unsigned aux = 0;
				for (int k = 0; k < size; ++k)
				{
					if(k != j && k != i)
					{
						++aux;
						if(type == 0 || type == 3)
						{
							outputFile << "G" << l1 << "_" << c1;
							outputFile << "c" << i << "_" << k;	
						}
						else
						{
							outputFile << "G" << l1 << "_" << c1;
							outputFile << "c" << k << "_" << i;
						}
						if(aux == size-2)
							outputFile << ";";
						else	
							outputFile << " + ";	
					}
				}
				outputFile << std::endl;	
			}
			else if(i != j && type == 4)
			{
				outputFile << "G" << l0 << "_" << c0;
				outputFile << "c" << i << "_" << j;
				outputFile << " <= ";
				outputFile << "G" << l1 << "_" << c1;
				outputFile << "c" << j << "_" << i;
				outputFile << ";" << std::endl;
			}
		}
	}
}

bool empty_line(matrix& cnot, unsigned line)
{
	for (int i = 0; i < cnot.size(); ++i)
	{
		if( cnot[i][line] > 0 )
			return false;
		if( cnot[line][i] > 0 )
			return false;
	}
	return true;
}

void writeDepSingle( matrix& output, unsigned l, unsigned c, unsigned size, unsigned type )
{
	if(type == 5)
		outputFile << "/* Writing Single Control dependency (" << l << "," << c << ") */" << std::endl;
	else if(type == 6)
		outputFile << "/* Writing Single Target dependency (" << l << "," << c  << ") */" << std::endl;
	else
		outputFile << "/* ERRORRRRRRRRRRRRRRRRR */" << std::endl;

	for (int i = 0; i < size; ++i)
	{
		// if( i == l || i == c)
		// {
			outputFile << "1";
			for (int j = 0; j < size; ++j)
			{
				if(i != j)
					outputFile << " - G" << l << "_" << c;
				if(i != j && type == 5)
					outputFile << "c" << i << "_" << j;
				else if(i != j && type == 6)
					outputFile << "c" << j << "_" << i;
			}
			outputFile << " <= ";
			for (int m = 0; m < size; ++m)
			{
				for (int n = 0; n < size; ++n)
				{
					if(m == l && n == c)
					{}	//nothing;
					else if(output[m][n] > 0)
					{
						for (int j = 0; j < size; ++j)
						{
							if(i != j)
							{
								outputFile << "G" << m << "_" << n;
								if(type == 5)
									outputFile << "c" << i << "_" << j << " + ";
								else if(i != j && type == 6)
									outputFile << "c" << j << "_" << i << " + ";

								outputFile << "G" << m << "_" << n;
								if(type == 5)
									outputFile << "c" << j << "_" << i;
								else if(i != j && type == 6)
									outputFile << "c" << i << "_" << j;

								outputFile << " + ";
							}
						}
					}
				}
			}
			outputFile << "0;" << std::endl;
		// }
	}
}

// Create the ghost gates
void ghostConnections(matrix& cnot)
{
	unsigned single;
	std::vector<int> ghost, a, b;
	// Search for single control or target
	for (int i = 0; i < cnot.size(); ++i)
	{
		single = 0;
		for (int j = 0; j < cnot.size(); ++j)
		{
			if( cnot[i][j] > 0 )
				++single;
			if( cnot[j][i] > 0)
				++single;
		}
		if( single == 1)
			ghost.push_back(i);
	}

	// for (int i = 0; i < ghost.size(); ++i)
	// {
	// 	std::cout << " " << ghost[i];
	// }
	// std::cout << std::endl;

	for (int i = 0; i < ghost.size(); ++i)
	{
		for (int j = i + 1; j < ghost.size(); ++j)
		{
			if( ghost[i] >= 0 && ghost[j] >= 0 && cnot[ghost[i]][ghost[j]] == 0 && cnot[ghost[j]][ghost[i]] == 0 )
			{
				for (int k = 0; k < cnot.size(); ++k)
				{
					if( cnot[ghost[i]][k] > 0 )
					{
						// std::cout << "Create ghost gate: " << ghost[i] << "-" << ghost[j] << " " << ghost[i] << "-" << k << std::endl;
						writeDep(ghost[i], k, ghost[i], ghost[j], cnot.size(), 0);
						a.push_back(ghost[i]);
						b.push_back(ghost[j]);
					}
					else if( cnot[k][ghost[i]] > 0 )
					{
						// std::cout << "Create ghost gate: " << ghost[i] << "-" << ghost[j] << " " << k << "-" << ghost[i] << std::endl;
						writeDep(k, ghost[i], ghost[i], ghost[j], cnot.size(), 3);
						a.push_back(ghost[i]);
						b.push_back(ghost[j]);
					}
					if( cnot[ghost[j]][k] > 0 )
					{
						// std::cout << "Create ghost gate: " << ghost[i] << "-" << ghost[j] << " " << ghost[j] << "-" << ghost[k] << std::endl;
						writeDep(ghost[j], k, ghost[i], ghost[j], cnot.size(), 2);
						a.push_back(ghost[i]);
						b.push_back(ghost[j]);
					}
					else if( cnot[k][ghost[j]] > 0 )
					{
						// std::cout << "Create ghost gate: " << ghost[i] << "-" << ghost[j] << " " << k << "-" << ghost[j] << std::endl;
						writeDep(k, ghost[j], ghost[i], ghost[j], cnot.size(), 1);
						a.push_back(ghost[i]);
						b.push_back(ghost[j]);
					}		
				}
				ghost[j] = -1;
				ghost[i] = -1;
				break;
			}
		}
	}
	for (int i = 0; i < a.size(); ++i)
		cnot[a[i]][b[i]] = 1;
}

// Naive solution
void getAllCombinations(matrix& output)
{
	// std::cout << "Getting all the dependencies..." << std::endl;
	enum type { cc, tt, ct, tc, in, sc, st };
	for (int i = 0; i < output.size(); ++i)
	{
		for (int j = 0; j < output[i].size(); ++j)
		{
			for (int m = 0; m < output.size(); ++m)
			{
				for (int n = 0; n < output[m].size(); ++n)
				{
					if(i != j && m != n && output[i][j] > 0 && output[m][n] > 0)
					{
						if(i == m && j != n)
							writeDep( i, j, m, n, output.size(),  cc);
						else if(i != m && j == n)
							writeDep( i, j, m, n, output.size(),  tt);
						else if(i == n && j != m)
							writeDep( i, j, m, n, output.size(),  ct);
						else if(i != n && j == m)
							writeDep( i, j, m, n, output.size(),  tc);
						else if(i == n && j == m)
							writeDep( i, j, m, n, output.size(),  in);
					}
				}
			}		
		}
	}
	// std::cout << "Done!" << std::endl;
}

bool alex_command::execute()
{
	circuit circ = env->store<circuit>().current();
	matrix output;
	// matrix qx4 = {{0,4,10,20,19,29,39,51,61,64,54,42,30,20,10,4},
	// 				{0,0,0,3,9,19,29,41,51,61,53,41,29,19,9,10},
	// 				{10,4,0,0,3,22,32,44,54,64,54,42,30,20,3,4},
	// 				{20,14,4,0,0,10,20,32,42,52,44,32,20,10,0,10},
	// 				{30,24,14,4,0,4,14,20,30,40,32,20,14,4,10,20},
	// 				{42,30,20,10,0,0,4,10,20,30,22,10,4,10,22,30},
	// 				{54,42,32,22,3,0,0,0,10,20,3,0,10,22,34,42},
	// 				{64,52,42,32,22,10,4,0,4,10,0,10,20,32,44,52},
	// 				{76,64,54,44,34,22,10,0,0,4,3,20,30,42,54,64},
	// 				{66,74,64,54,44,32,20,3,0,0,0,10,20,32,44,54},
	// 				{54,62,52,42,32,20,14,4,10,4,0,4,14,20,32,42},
	// 				{44,52,42,32,22,10,4,10,20,10,0,0,4,10,22,32},
	// 				{34,42,32,22,3,0,10,22,32,22,3,0,0,0,3,22},
	// 				{22,30,20,10,0,10,20,32,42,32,22,10,4,0,0,10},
	// 				{10,20,10,4,10,20,30,42,52,42,32,20,14,4,0,4},
	// 				{0,10,0,3,9,19,29,41,51,54,44,32,20,10,0,0}};
	// matrix qx4 ={{0,4,4,10,10},{0,0,4,10,10},{0,0,0,4,0},{3,3,0,0,0},{10,10,4,4,0}};
	matrix qx4 ={	{0,4,4,10,10},
					{0,0,4,10,10},
					{0,0,0,4,0},
					{12,12,0,0,0},
					{10,10,4,4,0}};

	if(filename.empty())
	{
		std::cout << "Missing output file. Use -f" << std::endl;
		return true;
	}
	
	outputFile.open (filename);

  	createMatrix( output, circ.lines() );
  	generateMatrixCnots( circ, output );
	// printMatrixCnots( output );
	printObjectiveFunction( qx4, output, getNumberDifGates(output) );
	if(getNumberDifGates(output) > 1)
		getAllCombinations(output);
	ghostConnections(output);
	printFirstRestriction( output );
	printEndRestriction( qx4, output, getNumberDifGates( output ) );
	printIntegerVariables( qx4, output, getNumberDifGates( output ) );
  	outputFile.close();
  	filename.clear();
	return true;
}

command::log_opt_t alex_command::log() const
{
  return log_opt_t({
	  {"runtime",       statistics->get<double>( "runtime" )}
	});
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
