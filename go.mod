// The oblikovati-traceon add-in: a c-shared library (.so/.dll) loaded by the host
// at runtime, integrating Traceon-equivalent electron-optics simulation (radially
// symmetric Boundary Element Method electrostatic/magnetostatic solver + charged
// particle tracer). The numerical core (./core) is a pure-Go, oracle-verified port
// of the upstream MPL-2.0 Traceon library and has NO host dependency; the engine
// (./engine) pulls geometry/materials from the host over the Apache-2.0 API, runs
// the core solver+tracer, and renders fields/trajectories back as client graphics.
//
// The SHIPPED library links only the Apache-2.0 contract (oblikovati.org/api) plus
// pure-Go numerics (gonum). The runtime boundary to the host is the C ABI, not Go
// (see ./include/oblikovati_addin.h). Sibling repos are resolved by the go.work at
// this repo's root (no committed replace); CI injects the equivalent replaces.
module oblikovati.org/traceon

go 1.24.0

// gonum.org/v1/gonum (pure-Go LU/QR/splines, wrapped behind core/linalg) is added in
// the PBI that first imports it; declaring it before use makes `go mod tidy` strip it.
require oblikovati.org/api v0.91.0
