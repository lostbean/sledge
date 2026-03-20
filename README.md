# sledge

Crystallographic texture analysis for polycrystalline materials.

## What is this?

When a metal is processed (rolled, forged, heat-treated), its crystal grains develop preferred orientations -- this is called *texture*. Texture profoundly affects material properties: strength, ductility, magnetic behavior, and corrosion resistance all depend on how grains are oriented. `sledge` is a Haskell library for representing, analyzing, and modeling crystallographic texture.

## What does it provide?

### Orientation representation

A crystal orientation is a rotation that maps the crystal's coordinate system onto the sample's coordinate system. Different communities and software tools prefer different representations. `sledge` supports all five standard representations with full interconversion:

- **Quaternions** -- the primary internal representation. Unit quaternions (4D vectors of length 1) offer singularity-free rotation composition and smooth interpolation (SLERP). They are compact and numerically well-behaved.
- **Euler angles** (Bunge convention) -- three angles (phi1, PHI, phi2) describing sequential rotations around Z, X', Z''. The most common representation in metallurgical literature, but prone to gimbal lock.
- **Axis-angle** -- a rotation axis (unit vector) and an angle. Intuitive but not unique.
- **Rodrigues vectors** -- the rotation axis scaled by the tangent of half the angle. Useful because the fundamental zone (the region of unique orientations) is a polyhedron in Rodrigues space.
- **Rotation matrices** -- 3x3 orthonormal matrices. Straightforward for vector rotation, but redundant (9 numbers for 3 degrees of freedom).

### Crystal symmetry

Crystal symmetry means that many different rotations produce the same physical crystal orientation. For example, cubic crystals (like iron, aluminum, copper) have 24 symmetry operators -- meaning each physical orientation corresponds to 24 equivalent quaternions. `sledge` provides symmetry group definitions (cubic, hexagonal, custom), fundamental zone computation, and misorientation calculation accounting for symmetry.

### Orientation Distribution Functions (ODFs)

An ODF describes the probability density of crystal orientations in a polycrystal. Given a set of measured orientations (e.g., from an EBSD scan of thousands of grains), `sledge` builds an ODF by placing Gaussian kernels at each measured orientation on a discrete grid in Rodrigues-Frank space. A vantage-point tree (from `queryforest`) indexes the grid for efficient neighbor lookup, so only nearby points (within 3 standard deviations) contribute to each kernel.

### Bingham distributions

The Bingham distribution is a probability distribution on the 3-sphere (the space of unit quaternions) that naturally models crystallographic texture. It is parameterized by three concentration values and three direction vectors, which together describe an ellipsoidal probability density on the orientation sphere.

- **Fitting** via maximum likelihood estimation: compute the scatter matrix of observed quaternions, eigendecompose to get principal directions, then optimize concentration parameters using BFGS
- **Normalization** via the saddle-point approximation of Kume and Wood (2005/2007) -- the Bingham normalization constant has no closed form, so a numerical approximation is necessary
- **Sampling** via Metropolis-Hastings with an angular central Gaussian proposal distribution

### Inverse Pole Figure (IPF) coloring

Assigns RGB colors to crystal orientations based on which crystallographic direction aligns with a chosen sample direction. This is the standard coloring scheme used in EBSD maps -- e.g., in a cubic IPF-Z map, grains with [001] parallel to the sample normal appear red, [101] green, and [111] blue.

### EBSD file I/O

Reads and writes the two main EBSD data file formats:

- **ANG** (TSL/OIM by EDAX) -- Euler angles in radians
- **CTF** (Channel Text File by Oxford/HKL) -- Euler angles in degrees

Bidirectional conversion between ANG and CTF is supported. Auto-detection reads either format.

### Sphere grids and projections

- **Tesseract grid** -- discretization of quaternion space using a 4-cube decomposition for binning orientation data
- **Icosahedral geodesic grid** -- uniform discretization of the 2-sphere for directional data, built by recursive subdivision of an icosahedron
- **Sphere projections** -- Lambert equal-area and stereographic projections for 2D visualization of pole figures

## Example

```haskell
import Texture.Orientation
import Texture.Bingham
import File.EBSD

-- Read EBSD data and fit a Bingham distribution
ebsd <- readEBSD "sample.ang"
let qs = case ebsd of
      Left ang  -> V.map rotation (nodes ang)
      Right ctf -> V.map rotation (nodes ctf)

let dist = fitBingham (U.convert qs)           -- fit Bingham to observed orientations
    density = binghamPDF dist someQuaternion    -- evaluate probability density
samples <- sampleBingham dist 10000             -- draw 10,000 samples
```

## Where is it used?

- **VirMat** -- the virtual microstructure generator uses `sledge` to assign crystallographic orientations to generated grains (via Bingham distribution sampling), compute IPF colors for visualization, and export ANG files for analysis in EBSD software (OIM, MTEX, etc.)

## How to build

```bash
# With Nix (recommended)
nix develop
cabal build --allow-newer

# With Cabal
cabal build

# Run tests
cabal test
```

## References

- A. Kume and A.T.A. Wood, "Saddlepoint approximations for the Bingham and Fisher-Bingham normalising constants," Biometrika, 2005.
- A. Kume and A.T.A. Wood, "On the derivatives of the normalising constant of the Bingham distribution," Statistics & Probability Letters, 2007.
- K. Shoemake, "Uniform random rotations," Graphics Gems III, 1992.
- F.L. Markley, "Averaging Quaternions," 2007.

## License

MIT -- see [LICENSE](./LICENSE).
