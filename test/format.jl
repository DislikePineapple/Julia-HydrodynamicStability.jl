using Test
using JuliaFormatter

# Format all files in the directory
@test format_file(".")
