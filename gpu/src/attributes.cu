/*
** copyright banner
*/

#include <fstream>

#include "attributes.h"
#include "check-return.h"

#define CSV_IO_NO_THREAD
#include "csv.h"    // https://github.com/ben-strasser/fast-cpp-csv-parser

void load_attribute( const std::string& attribute_name, const std::string& attribute_file, SimulationState& state );

/****
**
** Main function, load attributes as specified in a key:value pair text file.
**
****/

void load_attributes(const std::string& filename, SimulationState& state)
{
    if ( filename.size() == 0 )
    {
        printf( "No attributes specified.\n" );
        return;
    }

    // separate columns with a ':' (looks like key:value pairs), allow comment lines prefixed with '#'
    io::CSVReader<2, io::trim_chars<' ', '\t'>, io::no_quote_escape<':'>, io::throw_on_overflow, io::single_line_comment<'#'>> reader(filename);
    // attribute files do not have a header, the schema is these two columns:
    reader.set_header( "attribute", "attribute_file" );
    std::string attribute, attribute_file;

    while( reader.read_row( attribute, attribute_file ) )
    {
        load_attribute( attribute, attribute_file, state );
    }
}

// Forward declare helper functions which handle error checking.
void open_file( std::ifstream& file, const std::string& filename );
void read_data( void* buffer, size_t element_size, size_t element_count, std::ifstream& file, const std::string& filename );

/****
**
** Primary helper function for loading each attribute file.
**
****/

void load_attribute( const std::string& attribute_name, const std::string& filename, SimulationState& state )
{
    printf( "Loading %d '%s' attributes from file '%s'.\n", state.population_size, attribute_name.c_str(), filename.c_str() );
    std::ifstream file;
    open_file( file, filename );

    ATTR_TYPE* memory(nullptr);
    CheckReturn(cudaMallocManaged( &memory, state.population_size * sizeof(ATTR_TYPE) ));

    read_data( memory, sizeof(ATTR_TYPE), state.population_size, file, filename );
    state.attributes[attribute_name] = memory;
}

// Helper functions which handle error checking.

void open_file( std::ifstream& file, const std::string& filename )
{
    file.open( filename, std::ios::in | std::ios::binary );
    if ( !file )
    {
        fprintf( stderr, "Could not open '%s' for reading.\n", filename.c_str());
        exit( 1 );
    }
}

void read_data( void* buffer, size_t element_size, size_t element_count, std::ifstream& file, const std::string& filename )
{
    size_t byte_count = element_size * element_count;
    file.read( (char*)buffer, byte_count );
    if ( !file )
    {
        fprintf( stderr, "Failed to read %ld bytes (%ld elements) from attribute file '%s'.\n", byte_count, byte_count / element_size, filename.c_str() );
        exit( 1 );
    }
}
