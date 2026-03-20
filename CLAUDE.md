# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

sledge is a Haskell library providing basic tools for texture analysis in polycrystal materials — rotations, Bingham distributions, spherical harmonics, and EBSD file I/O (ANG, CTF formats).

## Build & Development

**Prerequisites:** Nix (provides GHC, Cabal, HLS, formatters).

```bash
nix develop                    # Enter dev shell
cabal build --allow-newer      # Build the library (deps fetched from GitHub via cabal.project)
```

## CI & Quality Checks

```bash
bash scripts/ci.sh             # Full CI: format check + build with -Wall -Werror + tests
bash scripts/pre-commit.sh     # Apply formatting to all files
```

CI runs format check, compilation with `-Wall -Werror` (fail on warnings), and tests.

## Testing

```bash
cabal test --allow-newer       # Run test suite (Tasty + QuickCheck + HUnit)
```

## Formatting

Three formatters via treefmt/nix:
- **Fourmolu** — Haskell (.hs)
- **cabal-fmt** — Cabal files (.cabal)
- **nixpkgs-fmt** — Nix files (.nix)

```bash
nix fmt                    # Format all files
nix fmt -- --ci            # Check formatting without modifying (CI mode)
```

Pre-commit hooks via Lefthook auto-format staged files.

## Git Conventions

- **No footer on commit messages.** Do not add `Co-Authored-By` or any other footer lines.
- Keep commit messages concise and descriptive.
