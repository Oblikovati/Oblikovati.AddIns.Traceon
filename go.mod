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

require (
	gonum.org/v1/gonum v0.17.0 // pure-Go dense linalg, wrapped behind core/linalg
	oblikovati.org/api v0.97.0
)
