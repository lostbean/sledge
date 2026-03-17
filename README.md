# sledge

**Basic tools for crystallographic texture analysis in polycrystalline materials**

sledge is a Haskell library for metallurgical texture research. It provides orientation representation and operations, texture representation via kernel density estimation and Bingham distributions, crystal symmetry handling, EBSD data file I/O, and sampling algorithms for orientation distribution functions (ODFs).

## Overview

sledge addresses three core problems in computational texture analysis:

1. **Orientation representation** -- converting between quaternions, Euler angles, axis-angle pairs, Rodrigues vectors, and rotation matrices, with full support for composition, inversion, misorientation, and interpolation.
2. **Texture modeling** -- constructing orientation distribution functions (ODFs) and direction distribution functions (DDFs) from EBSD measurements using kernel density estimation on vantage-point trees, or fitting parametric Bingham distributions via maximum likelihood estimation.
3. **EBSD data I/O** -- reading, writing, and converting between the ANG (TSL/OIM) and CTF (Oxford/HKL) file formats, as well as reading discrete ODF grid files.

## Modules

### Orientation and Rotation

#### `Texture.Orientation`

The central module. Defines orientation representations and operations following the convention e1 -> RD, e2 -> TD, e3 -> ND.

**Key types:**

- `Quaternion` -- Unit quaternion (the primary internal representation). Unboxed-vector derivable, with a `Random` instance using the Shoemake uniform random rotation algorithm.
- `Euler` -- Bunge-convention Euler angles (phi1, PHI, phi2) with sequence Z-X'-Z''.
- `AxisPair` -- Axis-angle representation (normalized direction + angle in radians).
- `Rodrigues` -- Frank-Rodrigues vector representation.
- `RotMatrix` -- 3x3 orthonormal rotation matrix.
- `Deg`, `Rad` -- Type-safe angle wrappers with an `Angle` typeclass for conversion.

**Key typeclass:**

- `Rot a` -- Defines `(#<=)` (composition), `(-#-)` (crystal-frame misorientation), `(-@-)` (sample-frame misorientation), `invert`, `zerorot`, `toQuaternion`, `fromQuaternion`, and `getOmega`. All five representations are instances.

**Key functions:**

- `activeVecRotation`, `passiveVecRotation` -- Rotate a vector actively or passively by a quaternion.
- `quaterInterpolation` -- SLERP interpolation between two quaternions.
- `averageQuaternion` -- Mean orientation using Markley's eigenvector method on the scatter matrix.
- `composeQ0` -- Fast computation of only the q0 component (used for fast misorientation angle queries).

#### `Texture.Symmetry`

Crystal symmetry operations and fundamental zone computation.

**Key types:**

- `Symm` -- Crystal symmetry group: `Cubic`, `Hexagonal`, or `Custom`.
- `SymmOp` -- A symmetry operator (wraps a `Quaternion`).

**Key functions:**

- `getSymmOps` -- Generate all symmetry operators for a symmetry group (e.g., 24 for cubic).
- `toFZ` -- Map a quaternion into the fundamental zone (minimum rotation angle among symmetric equivalents).
- `getMisoAngle` -- Minimum misorientation angle between two orientations accounting for crystal symmetry.
- `getAllSymmVec`, `getUniqueSymmVecs` -- All (or unique) symmetric equivalents of a direction vector.
- `isInRodriFZ` -- Fundamental zone test using Rodrigues-Frank space planes (faster for grid filtering).

### Texture Distributions

#### `Texture.ODF`

Orientation Distribution Function built by kernel density estimation on a discrete grid in Rodrigues-Frank fundamental zone space, backed by a vantage-point tree for fast neighbor queries.

**Key type:**

- `ODF` -- Contains the grid of quaternions, intensity values, a VP-tree index, symmetry, kernel width, grid size, and step.

**Key functions:**

- `buildEmptyODF` -- Construct an empty ODF given kernel width, symmetry, and grid step size.
- `addPoints` -- Accumulate orientation data using Gaussian kernel density estimation.
- `getODFeval` -- Return a function that evaluates the ODF intensity at any quaternion.
- `getMaxOrientation` -- Find the orientation with maximum intensity.
- `renderODFVTK` -- Export the ODF to VTK format for visualization.

#### `Texture.DDF`

Direction Distribution Function -- kernel density estimation on the unit sphere (S2) for directional data (e.g., pole figures). Uses an icosahedral geodesic grid (`IsoSphere`) and a VP-tree.

#### `Texture.Kernel`

Gaussian kernel density estimation engine for both ODF and DDF. Uses mutable ST vectors internally and the VP-tree for efficient neighbor lookup (only points within 3 sigma contribute).

#### `Texture.Bingham`

The Bingham distribution on the 3-sphere (S3), used for parametric modeling of crystallographic texture.

**Key type:**

- `Bingham` -- Contains concentration values, direction matrix, normalization constant F, mode quaternion, partial derivatives, and scatter matrix.

**Key functions:**

- `mkBingham` -- Construct from three (concentration, direction) pairs.
- `binghamPDF` -- Evaluate the probability density function at a quaternion.
- `fitBingham` -- Fit a Bingham distribution to observed quaternions using MLE (scatter matrix eigendecomposition + BFGS optimization).
- `sampleBingham` -- Generate samples using Metropolis-Hastings with an angular central Gaussian proposal.
- `bingProduct` -- Product of two Bingham distributions.

#### `Texture.Bingham.Constant` (internal)

Computation of the Bingham normalization constant via saddle-point approximation (Kume and Wood, 2005/2007).

### Sampling

#### `Texture.Sampler`

Hit-and-Run Slice Sampler -- a general-purpose MCMC sampler for arbitrary probability distributions, applicable to quaternions, 2D vectors, and scalars. Features adaptive step sizes and long-range exploration shots.

### Grids and Geometry

#### `Texture.TesseractGrid`

Discretization of quaternion space using a tesseract (4-cube) decomposition covering S3 with four cubic charts.

- `genQuaternionGrid` -- Generate a uniform grid of quaternions.
- `binningTesseract`, `accuTesseract` -- Bin quaternion data onto the grid.
- `plotTesseract` -- VTK export.

#### `Texture.IsoSphere`

Geodesic grid on the 2-sphere based on recursive subdivision of an icosahedron.

- `isoSphere` -- Create at a given subdivision level.
- `nearestPoint` -- Fast hierarchical nearest-point query.
- `genIsoSphereSO3Grid` -- Generate SO3 grid using icosahedral shells filtered to the Rodrigues FZ.

#### `Texture.HyperSphere`

Coordinate systems and grids on S2 (SO2) and S3 (SO3).

- `so3ToQuaternion`, `quaternionToSO3` -- Convert between SO3 coordinates and quaternions.
- `getSO2Grid`, `getSO3Grid` -- Regular grids on spheres.
- VTK rendering for scalar fields on spheres.

#### `Texture.SphereProjection`

2D projections: Lambert equal-area and stereographic.

### Inverse Pole Figure

#### `Texture.IPF`

IPF color computation. Assigns RGB colors to orientations based on the crystallographic direction parallel to a chosen sample direction.

- `getIPFColor` -- Compute IPF color for a quaternion in a given symmetry group.
- `genIPFLegend` -- Generate an IPF color legend.

### File I/O

#### `File.ANGReader` / `File.ANGWriter`

Read and write ANG files (TSL/OIM EBSD format). Euler angles parsed in radians, immediately converted to quaternions.

- `parseANG` -- Read and parse an ANG file.
- `renderANGFile` -- Write an ANG file.
- `angToVoxBox` -- Convert ANG data to a `VoxBox`.

#### `File.CTFReader` / `File.CTFWriter`

Read and write CTF files (Oxford/HKL Channel Text File format). Euler angles parsed in degrees.

- `parseCTF` -- Read and parse a CTF file.
- `renderCTFFile` -- Write a CTF file.
- `ctfToVoxBox` -- Convert CTF data to a `VoxBox`.

#### `File.EBSD`

Unified EBSD interface with format conversion and auto-detection.

- `readEBSD` -- Auto-detect and parse either ANG or CTF format.
- `readEBSDToVoxBox` -- Read EBSD data and convert to a `VoxBox`.
- `writeANG`, `writeCTF` -- Write any EBSD data to ANG or CTF format.
- `EBSD` typeclass with `toANG` and `toCTF` for bidirectional conversion.

#### `File.ODFReader`

Parser for discrete ODF grid files (text format).

## File Formats

| Format | Extension | Read | Write | Angles | Origin |
|--------|-----------|------|-------|--------|--------|
| ANG    | `.ang`    | Yes  | Yes   | Radians (Euler) | TSL/OIM (EDAX) |
| CTF    | `.ctf`    | Yes  | Yes   | Degrees (Euler) | Oxford/HKL (Channel) |
| ODF    | (text)    | Yes  | No    | Discrete grid   | Generic |

Bidirectional ANG <-> CTF conversion is supported with approximate field mapping.

## Algorithms

### Kernel Density Estimation

The ODF and DDF are built by placing Gaussian kernels at each observed orientation/direction. A vantage-point tree indexes the grid points for efficient neighbor lookup (only points within 3 sigma contribute).

### Bingham Distribution Fitting (MLE)

1. Compute the scatter matrix of the observed quaternions.
2. Eigendecompose to obtain principal directions and eigenvalues.
3. Find concentration parameters via BFGS optimization matching theoretical and observed partial derivatives of the normalization constant.

### Bingham Normalization (Saddle-Point Approximation)

Uses the method of Kume and Wood (2005, 2007): find the saddle point via Newton's method, apply the Laplace-type expansion with a fourth-order correction.

### Metropolis-Hastings Sampling (Bingham)

Proposal distribution: angular central Gaussian parameterized by the Bingham scatter matrix eigenvectors. Standard MH acceptance ratio. Burn-in of 20 samples.

### Hit-and-Run Slice Sampling

General MCMC sampler: draw a vertical level below current density, choose a random direction, expand an interval until both endpoints fall below the level, sample uniformly with shrinking on rejection. Adaptive step sizes and periodic long-range exploration shots.

### Fundamental Zone Computation

Two approaches: quaternion-based (apply all symmetry operators, maximize |q0|) and Rodrigues-Frank plane test (faster for bulk grid filtering).

## Usage Examples

### Orientation manipulation

```haskell
import Texture.Orientation

let e1 = mkEuler (Deg 30) (Deg 45) (Deg 60)
    q1 = toQuaternion e1

-- Compose rotations
let q3 = q1 #<= q2

-- Misorientation in crystal frame
let miso = q2 -#- q1

-- Convert between representations
let rod = fromQuaternion q1 :: Rodrigues
    mat = fromQuaternion q1 :: RotMatrix
```

### Building an ODF from EBSD data

```haskell
import File.EBSD
import Texture.ODF

ebsd <- readEBSD "sample.ang"
let qs = case ebsd of
      Left ang  -> V.map rotation (nodes ang)
      Right ctf -> V.map rotation (nodes ctf)

let odf = addPoints (U.convert qs) (buildEmptyODF (Deg 2.5) Cubic (Deg 3))
```

### Fitting a Bingham distribution

```haskell
import Texture.Bingham

let dist = fitBingham (U.fromList quaternionList)
let density = binghamPDF dist someQuaternion
samples <- sampleBingham dist 10000
```

### File format conversion

```haskell
import File.EBSD

ebsd <- readEBSD "input.ang"
case ebsd of
  Left ang -> writeCTF "output.ctf" ang
  Right ctf -> writeCTF "output.ctf" ctf
```

## Dependencies

| Package | Purpose |
|---------|---------|
| `base` (>= 4, < 5) | Standard library |
| `attoparsec` | Fast parsing for ANG/CTF/ODF files |
| `binary` | Binary serialization (TesseractGrid) |
| `blaze-builder` | Efficient ByteString construction for file writers |
| `bytestring` | Byte string I/O |
| `deepseq` | NFData instances |
| `hammer` | VTK output, VoxBox, BFGS optimizer |
| `linear-vect` | Vector/matrix algebra, eigendecomposition |
| `primitive` | Mutable vector primitives |
| `queryforest` | Vantage-point tree for spatial queries |
| `random` | Random number generation |
| `text` | Text processing |
| `vector` | Boxed and unboxed vectors |
| `vector-binary-instances` | Binary instances for vectors |
| `vector-th-unbox` | Template Haskell for Unbox instances |

## Building

```bash
# With Stack (from the parent project)
stack build sledge

# With Cabal
cabal build

# With Nix
nix develop
cabal build

# Run tests
stack test sledge
```

Optional visual test executables (behind the `tests` flag):

```bash
stack build sledge --flag sledge:tests
# test-kernel-sampling, test-odf, test-sampler
```

## References

- Conversion of EBSD data by a quaternion based algorithm to be used for grain structure simulations.
- Orientation Library Manual by Romain Quey.
- K. Shoemake, "Uniform random rotations," in Graphics Gems III, 1992.
- F.L. Markley, "Averaging Quaternions," 2007.
- A. Kume and A.T.A. Wood, "Saddlepoint approximations for the Bingham and Fisher-Bingham normalising constants," Biometrika, 2005.
- A. Kume and A.T.A. Wood, "On the derivatives of the normalising constant of the Bingham distribution," Statistics & Probability Letters, 2007.

## Author

Edgar Gomes de Araujo (<talktoedgar@gmail.com>)

## License

MIT -- see [LICENSE](./LICENSE).
