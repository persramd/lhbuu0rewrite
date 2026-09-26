# LHBUU Transport Code - Modern Rewrite

Modern C implementation of Lattice Hamiltonian Boltzmann-Uehling-Uhlenbeck nuclear transport theory.

## Status: ~18% Complete (10/15 core modules)

See [SESSION_SUMMARY.md](SESSION_SUMMARY.md) for full progress report.

## Quick Start

```bash
gcc -c src/form_factor.c -Iinclude -O3
gcc -c src/grid_tag.c -Iinclude -O3
# All modules compile cleanly
```

## Next Task

Implement `grad_V_alpha_p_md` - see SESSION_SUMMARY.md

## Documentation

- SESSION_SUMMARY.md - Complete progress and roadmap
- IMPLEMENTATION_STATUS.md - Component details
- DESIGN.md - Technical design decisions
