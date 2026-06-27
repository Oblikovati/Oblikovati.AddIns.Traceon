// SPDX-License-Identifier: MPL-2.0

// Package elliptic provides the complete elliptic integrals K and E used by the
// radially-symmetric ring kernels. It is a faithful port of traceon/backend/elliptic.c,
// which uses the Chebyshev approximations of W. J. Cody (1965) augmented with the
// argument-reduction tricks from the SciPy ellipkm1 documentation.
//
// Convention (matching the upstream and SciPy): the argument is the parameter m = k^2,
// NOT the modulus k. K(m) = integral_0^{pi/2} dθ / sqrt(1 - m sin^2 θ).
package elliptic

import "math"

// cody1965 evaluates the degree-7 Chebyshev series sum_i (A[i] + L*B[i]) * p^i with
// L = log(1/p). Shared by the km1/em1 forms, which differ only in their coefficients.
func cody1965(p float64, a, b *[8]float64) float64 {
	l := math.Log(1.0 / p)
	sum := 0.0
	pi := 1.0 // p^i, built incrementally to avoid pow() per term (matches C loop order)
	for i := 0; i < 8; i++ {
		sum += (a[i] + l*b[i]) * pi
		pi *= p
	}
	return sum
}

// kCoeffA / kCoeffB are Cody's coefficients for K; A[0] = log(4).
var (
	kCoeffA = [8]float64{
		math.Log(4.0),
		9.65736020516771e-2, 3.08909633861795e-2, 1.52618320622534e-2,
		1.25565693543211e-2, 1.68695685967517e-2, 1.09423810688623e-2, 1.40704915496101e-3,
	}
	kCoeffB = [8]float64{
		1.0 / 2.0,
		1.24999998585309e-1, 7.03114105853296e-2, 4.87379510945218e-2,
		3.57218443007327e-2, 2.09857677336790e-2, 5.81807961871996e-3, 3.42805719229748e-4,
	}
	eCoeffA = [8]float64{
		1,
		4.43147193467733e-1, 5.68115681053803e-2, 2.21862206993846e-2,
		1.56847700239786e-2, 1.92284389022977e-2, 1.21819481486695e-2, 1.55618744745296e-3,
	}
	eCoeffB = [8]float64{
		0,
		2.49999998448655e-1, 9.37488062098189e-2, 5.84950297066166e-2,
		4.09074821593164e-2, 2.35091602564984e-2, 6.45682247315060e-3, 3.78886487349367e-4,
	}
)

// Ellipkm1 is the complete elliptic integral K evaluated at m = 1 - p, i.e. K(1-p),
// accurate for small p (near the m→1 singularity). Port of ellipkm1.
func Ellipkm1(p float64) float64 { return cody1965(p, &kCoeffA, &kCoeffB) }

// Ellipem1 is the complete elliptic integral E evaluated at m = 1 - p, i.e. E(1-p).
// Port of ellipem1.
func Ellipem1(p float64) float64 { return cody1965(p, &eCoeffA, &eCoeffB) }

// Ellipk is the complete elliptic integral K(m) for parameter m. For m > -1 it reduces
// to Ellipkm1(1-m); otherwise it uses the reciprocal-modulus transform. Port of ellipk.
func Ellipk(m float64) float64 {
	if m > -1 {
		return Ellipkm1(1 - m)
	}
	return Ellipkm1(1.0/(1-m)) / math.Sqrt(m)
}

// Ellipe is the complete elliptic integral E(m) for parameter m. For 0 <= m <= 1 it
// reduces to Ellipem1(1-m); otherwise it uses the modulus transform. Port of ellipe.
func Ellipe(m float64) float64 {
	if 0 <= m && m <= 1 {
		return Ellipem1(1 - m)
	}
	return Ellipem1(-1/(m-1.0)) * math.Sqrt(1-m)
}
