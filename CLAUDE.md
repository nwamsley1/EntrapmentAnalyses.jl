# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Julia package called "EntrapmentAnalyses" (v0.1.0) that appears to be in early development. The package is set up with the standard Julia package structure.

## Common Development Commands

### Julia Package Management
```bash
# Enter Julia REPL
julia

# Activate the project environment
] activate .

# Install dependencies
] instantiate

# Update dependencies
] update

# Add a new dependency
] add PackageName

# Remove a dependency
] rm PackageName
```

### Running Tests
```bash
# From Julia REPL with activated environment
] test

# Or from command line
julia --project=. -e 'using Pkg; Pkg.test()'

# Run tests with coverage
julia --project=. --code-coverage=user -e 'using Pkg; Pkg.test()'
```

### Development Workflow
```bash
# Start Julia with Revise for auto-reloading
julia --project=.

# In the REPL
using Revise
using EntrapmentAnalyses
```

## Code Architecture

The package follows standard Julia package conventions:

- `src/EntrapmentAnalyses.jl`: Main module file that exports the public API
- `test/runtests.jl`: Test suite entry point
- `Project.toml`: Package metadata and dependencies
- `Manifest.toml`: Locked dependency versions

## Important Notes

1. **Test File Issue**: The test file (`test/runtests.jl`) currently references `MyPackage` instead of `EntrapmentAnalyses`. This needs to be fixed before tests can run properly.

2. **Revise.jl Integration**: The package includes Revise as a dependency, which enables automatic code reloading during development. Always start Julia sessions with Revise loaded for a smoother development experience.

3. **Package UUID**: The package has UUID `60062a3e-e55a-4920-8913-bf16cdadf465` which uniquely identifies it in the Julia ecosystem.

## Julia-Specific Conventions

- Use `@test` macro for unit tests
- Follow Julia naming conventions: modules use PascalCase, functions use snake_case
- Export public API functions at the module level
- Include docstrings using triple quotes before function definitions
- Use type annotations for better performance and clarity when beneficial