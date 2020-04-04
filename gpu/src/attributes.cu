/*
** copyright banner
*/

#include <fstream>

#include "attributes.h"
#include "check-return.h"

#if !defined(NDEBUG)
static bool verbose = false;
#else
static bool verbose = true;
#endif

void check_file( const std::ifstream& file, const std::string& filename );
void split_line( const std::string& line, std::string& name, std::string& file );
void load_attribute( const std::string& attribute_name, const std::string& attribute_file, SimulationState& state );

void load_attributes(const std::string& filename, SimulationState& state)
{
    if ( filename.size() == 0 )
    {
        printf( "No attributes specified.\n" );
        return;
    }

    std::ifstream infile(filename);
    check_file( infile, filename );
    std::string line;
    while (std::getline(infile, line))
    {
        if ( verbose ) printf( "Read '%s' from '%s'.\n", line.c_str(), filename.c_str());
        std::string attribute_name, attribute_file;
        split_line( line, attribute_name, attribute_file );
        load_attribute( attribute_name, attribute_file, state );
    }
}

void check_file( const std::ifstream& file, const std::string& filename )
{
    if ( file.fail() )
    {
        fprintf( stderr, "Could not open '%s' for reading.\n", filename.c_str() );
        exit( 1 );
    }
}

void split_line( const std::string& line, std::string& name, std::string& file )
{
    size_t colon = line.find(':');
    if ( colon != std::string::npos )
    {
        name = line.substr( 0, colon );
        file = line.substr( colon + 1 );
        if ( verbose ) printf( "Extracted attribute '%s' and file '%s' from '%s'.\n", name.c_str(), file.c_str(), line.c_str() );
    }
    else
    {
        fprintf( stderr, "Attribute file should consist of name:file pairs, got '%s'.\n", line.c_str() );
        exit( 1 );
    }
}

FILE* open_file( const std::string& filename );
void read_data( void* buffer, size_t element_size, size_t element_count, FILE* file, const std::string& filename );

void load_attribute( const std::string& attribute_name, const std::string& attribute_file, SimulationState& state )
{
    printf( "Loading %d '%s' attributes from file '%s'.\n", state.population_size, attribute_name.c_str(), attribute_file.c_str() );
    auto fp = open_file( attribute_file );
    ATTR_TYPE* memory(nullptr);
    CheckReturn(cudaMallocManaged(&memory, state.population_size * sizeof(ATTR_TYPE)));
    read_data( memory, sizeof(ATTR_TYPE), state.population_size, fp, attribute_file );
    state.attributes[attribute_name] = memory;
}

FILE* open_file( const std::string& filename )
{
    auto fp = fopen( filename.c_str(), "r" );
    if ( !fp )
    {
        fprintf( stderr, "Could not open attribute file '%s' for reading.\n", filename.c_str() );
        exit( 1 );
    }

    return fp;
}

void read_data( void* buffer, size_t element_size, size_t element_count, FILE* file, const std::string& filename )
{
    size_t count = fread( buffer, element_size, element_count, file );
    if ( count != element_count )
    {
        fprintf( stderr, "Failed to read %ld elements from attribute file '%s' - got %ld.\n", element_count, filename.c_str(), count );
        exit( 1 );
    }
}
