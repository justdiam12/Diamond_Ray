/**
 * @file
 * @brief C++ main program for RAM
 */

#include "RamRunner.hh"

#include <iostream>
#include <string>
#include <vector>

/// Main program in C++ for Range-dependent Acoustic Model (RAM) 1.5

int main( int argc, char *argv[] )
{
    std::string inFile("ram.in");
    std::string lineFile("tl.line");
    std::string gridFile("tl.grid");
    
    // parse command line for input file 'ram.in'
    if (argc > 1)
    {
        size_t found;
        
        // if input param filename ends with '.in'
        //   replace suffix for output filenames
        // else just append suffixes.
        inFile = argv[1];
        found = inFile.rfind(".in");
        if (found!=std::string::npos)
        {
            lineFile = inFile;
            gridFile = inFile;
            lineFile.replace (found, 3, ".line");
            gridFile.replace (found, 3, ".grid");
        }
        else
        {
            lineFile = inFile + ".line";
            gridFile = inFile + ".grid";
        }
    }
    std::cout << "input file = "<< inFile << std::endl;
    std::cout << "TL vs range = " << lineFile << std::endl;
    std::cout << "TL field = " << gridFile << std::endl;
    
    // Create the object that will do all of the work
    RamRunner rar( inFile, lineFile, gridFile );
    
    /// Write the transmission loss files in the original RAM format
    rar.writeLossFiles();
    
    // March the acoustic field out in range.
    while ( rar.moreToCome() ) {
        
        // Propagate the field to the next range.
        (void)rar.PropagateToNextRange();

        /// Write the transmission loss files in the original RAM format
        rar.writeLossFiles();
        
    }   // end while rar.moreToCome()
    
    // Files are closed by RamRunner destructor
    return 0;
}

