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

go 1.27.0

require (
	gonum.org/v1/gonum v0.17.0 // pure-Go dense linalg, wrapped behind core/linalg
	oblikovati.org/api v0.154.0
)

require (
	github.com/bitfield/gotestdox v0.2.2 // indirect
	github.com/dnephin/pflag v1.0.7 // indirect
	github.com/fatih/color v1.18.0 // indirect
	github.com/fsnotify/fsnotify v1.8.0 // indirect
	github.com/google/shlex v0.0.0-20191202100458-e7afc7fbc510 // indirect
	github.com/mattn/go-colorable v0.1.13 // indirect
	github.com/mattn/go-isatty v0.0.20 // indirect
	golang.org/x/mod v0.25.0 // indirect
	golang.org/x/sync v0.15.0 // indirect
	golang.org/x/sys v0.33.0 // indirect
	golang.org/x/term v0.32.0 // indirect
	golang.org/x/text v0.23.0 // indirect
	golang.org/x/tools v0.34.0 // indirect
	gotest.tools/gotestsum v1.12.3 // indirect
)

tool gotest.tools/gotestsum
