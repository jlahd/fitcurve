// An Algorithm for Automatically Fitting Digitized Curves
// by Philip J. Schneider
// from "Graphics Gems", Academic Press, 1990
//
// Converted to golang from the original C

package fitcurve

import "math"

type (
	// Point2 is a 2D point.
	Point2 [2]float64

	// BezierCurve defines the ends and two control points of a cubic Bezier curve.
	BezierCurve [4]Point2
)

const maxIterations = 4

// FitCurve fits a Bezier curve to a set of digitized points.
func FitCurve(d []Point2, error float64) []BezierCurve {
	tHat1 := computeLeftTangent(d, 0)
	tHat2 := computeRightTangent(d, len(d)-1)
	return fitCubic(d, tHat1, tHat2, error)
}

// fitCubic fits a Bezier curve to a (sub)set of digitized points.
func fitCubic(d []Point2, tHat1, tHat2 Point2, error float64) []BezierCurve {
	iterationError := error * error

	// Use heuristic if region only has two points in it
	if len(d) == 2 {
		dist := v2DistanceBetween2Points(d[1], d[0]) / 3.0
		return []BezierCurve{
			{
				d[0],
				v2Add(d[0], v2Scale(tHat1, dist)),
				v2Add(d[1], v2Scale(tHat2, dist)),
				d[1],
			},
		}
	}

	// Parameterize points, and attempt to fit curve
	u := chordLengthParameterize(d)
	bezCurve := generateBezier(d, u, tHat1, tHat2)

	// Find max deviation of points to fitted curve
	maxError, splitPoint := computeMaxError(d, bezCurve, u)
	if maxError < error {
		return []BezierCurve{bezCurve}
	}

	// If error not too large, try some reparameterization
	// and iteration
	if maxError < iterationError {
		for i := 0; i < maxIterations; i++ {
			uPrime := reparameterize(d, u, bezCurve)
			bezCurve = generateBezier(d, uPrime, tHat1, tHat2)
			maxError, splitPoint = computeMaxError(d, bezCurve, uPrime)
			if maxError < error {
				return []BezierCurve{bezCurve}
			}
			u = uPrime
		}
	}

	// Fitting failed -- split at max error point and fit recursively
	tHatCenter := computeCenterTangent(d, splitPoint)
	c1 := fitCubic(d[:splitPoint+1], tHat1, tHatCenter, error)
	tHatCenter = v2Negate(tHatCenter)
	c2 := fitCubic(d[splitPoint:], tHatCenter, tHat2, error)
	return append(c1, c2...)
}

// generateBezier uses least-squares method to find Bezier control points for region.
func generateBezier(d []Point2, uPrime []float64, tHat1, tHat2 Point2) BezierCurve {
	last := len(d) - 1

	// Compute the A's
	a := make([][2]Point2, len(uPrime))
	for i, u := range uPrime {
		v1 := v2Scale(tHat1, b1(u))
		v2 := v2Scale(tHat2, b2(u))
		a[i][0] = v1
		a[i][1] = v2
	}

	// Create the C and X matrices
	var c [2][2]float64
	var x [2]float64

	for i, ai := range a {
		c[0][0] += v2Dot(ai[0], ai[0])
		c[0][1] += v2Dot(ai[0], ai[1])
		c[1][0] = c[0][1]
		c[1][1] += v2Dot(ai[1], ai[1])
		tmp := v2Sub(d[i],
			v2Add(
				v2Scale(d[0], b1(uPrime[i])),
				v2Add(
					v2Scale(d[0], b1(uPrime[i])),
					v2Add(
						v2Scale(d[last], b2(uPrime[i])),
						v2Scale(d[last], b3(uPrime[i]))))))
		x[0] += v2Dot(ai[0], tmp)
		x[1] += v2Dot(ai[1], tmp)
	}

	// Compute the determinants of C and X
	detC0C1 := c[0][0]*c[1][1] - c[1][0]*c[0][1]
	detC0X := c[0][0]*x[1] - c[1][0]*x[0]
	detXC1 := x[0]*c[1][1] - x[1]*c[0][1]

	// Finally, derive alpha values
	var alphaL, alphaR float64
	if detC0C1 != 0 {
		alphaL = detXC1 / detC0C1
		alphaR = detC0X / detC0C1
	}

	// If alpha negative, use the Wu/Barsky heuristic (see text)
	// (if alpha is 0, you get coincident control points that lead to
	// divide by zero in any subsequent newtonRaphsonRootFind() call.
	segLength := v2DistanceBetween2Points(d[last], d[0])
	epsilon := 1.0e-6 * segLength
	if alphaL < epsilon || alphaR < epsilon {
		// fall back on standard (probably inaccurate) formula, and subdivide further if needed.
		dist := segLength / 3.0
		return BezierCurve{
			d[0],
			v2Add(d[0], v2Scale(tHat1, dist)),
			v2Add(d[last], v2Scale(tHat2, dist)),
			d[last],
		}
	}

	// First and last control points of the Bezier curve are
	// positioned exactly at the first and last data points
	// Control points 1 and 2 are positioned an alpha distance out
	// on the tangent vectors, left and right, respectively
	return BezierCurve{
		d[0],
		v2Add(d[0], v2Scale(tHat1, alphaL)),
		v2Add(d[last], v2Scale(tHat1, alphaR)),
		d[last],
	}
}

// reparameterize tries, given set of points and their parameterization, to find
// a better parameterization.
func reparameterize(d []Point2, u []float64, bezCurve BezierCurve) []float64 {
	uPrime := make([]float64, len(d))
	for i, pt := range d {
		uPrime[i] = newtonRaphsonRootFind(bezCurve, pt, u[i])
	}
	return uPrime
}

// newtonRaphsonRootFind uses Newton-Raphson iteration to find better root.
func newtonRaphsonRootFind(q BezierCurve, p Point2, u float64) float64 {
	// Compute Q(u)
	qu := bezierII(3, q[:], u)

	// Generate control vertices for Q'
	var q1 [3]Point2
	for i := range q1 {
		q1[i][0] = (q[i+1][0] - q[i][0]) * 3.0
		q1[i][1] = (q[i+1][1] - q[i][1]) * 3.0
	}

	// Generate control vertices for Q''
	var q2 [2]Point2
	for i := range q2 {
		q2[i][0] = (q1[i+1][0] - q1[i][0]) * 2.0
		q2[i][1] = (q1[i+1][1] - q1[i][1]) * 2.0
	}

	// Compute Q'(u) and Q''(u)
	q1u := bezierII(2, q1[:], u)
	q2u := bezierII(1, q2[:], u)

	// Compute f(u)/f'(u)
	numerator := (qu[0]-p[0])*(q1u[0]) + (qu[1]-p[1])*(q1u[1])
	denominator := (q1u[0])*(q1u[0]) + (q1u[1])*(q1u[1]) +
		(qu[0]-p[0])*(q2u[0]) + (qu[1]-p[1])*(q2u[1])
	if denominator == 0.0 {
		return u
	}

	// u = u - f(u)/f'(u)
	return u - (numerator / denominator)
}

// bezier evaluates a Bezier curve at a particular parameter value
func bezierII(degree int, v []Point2, t float64) Point2 {
	// Copy array
	vtemp := make([]Point2, degree+1)
	copy(vtemp, v[:degree+1])

	// Triangle computation
	for i := 1; i <= degree; i++ {
		for j := 0; j <= degree-i; j++ {
			vtemp[j][0] = (1.0-t)*vtemp[j][0] + t*vtemp[j+1][0]
			vtemp[j][1] = (1.0-t)*vtemp[j][1] + t*vtemp[j+1][1]
		}
	}

	return vtemp[0]
}

// B0, B1, B2, B3 are Bezier multipliers
func b0(u float64) float64 {
	tmp := 1.0 - u
	return tmp * tmp * tmp
}

func b1(u float64) float64 {
	tmp := 1.0 - u
	return (3 * u * (tmp * tmp))
}

func b2(u float64) float64 {
	tmp := 1.0 - u
	return (3 * u * u * tmp)
}

func b3(u float64) float64 {
	return u * u * u
}

// computeLeftTangent, computeRightTangent, computeCenterTangent
// approximate unit tangents at endpoints and "center" of digitized curve
func computeLeftTangent(d []Point2, end int) Point2 {
	return v2Normalize(v2Sub(d[end+1], d[end]))
}

func computeRightTangent(d []Point2, end int) Point2 {
	return v2Normalize(v2Sub(d[end-1], d[end]))
}

func computeCenterTangent(d []Point2, center int) Point2 {
	v1 := v2Sub(d[center-1], d[center])
	v2 := v2Sub(d[center], d[center+1])
	return v2Normalize(Point2{
		(v1[0] + v2[0]) / 2,
		(v1[1] + v2[1]) / 2,
	})
}

// chordLengthParameterize assigns parameter values to digitized points
// using relative distances between points.
func chordLengthParameterize(d []Point2) []float64 {
	u := make([]float64, len(d))
	u[0] = 0.0
	for i := 1; i < len(d); i++ {
		u[i] = u[i-1] + v2DistanceBetween2Points(d[i], d[i-1])
	}

	for i := 1; i < len(d); i++ {
		u[i] = u[i] / u[len(u)-1]
	}

	return u
}

// computeMaxError finds the maximum squared distance of digitized points
// to fitted curve.
func computeMaxError(d []Point2, bezCurve BezierCurve, u []float64) (float64, int) {
	splitPoint := len(d) / 2
	maxDist := 0.0
	for i := 1; i < len(d)-1; i++ {
		p := bezierII(3, bezCurve[:], u[i])
		v := v2Sub(p, d[i])
		dist := v2SquaredLength(v)
		if dist >= maxDist {
			maxDist = dist
			splitPoint = i
		}
	}
	return maxDist, splitPoint
}

func v2DistanceBetween2Points(a, b Point2) float64 {
	dx := a[0] - b[0]
	dy := a[1] - b[1]
	return math.Sqrt(dx*dx + dy*dy)
}

func v2Add(a, b Point2) Point2 {
	return Point2{a[0] + b[0], a[1] + b[1]}
}

func v2Sub(a, b Point2) Point2 {
	return Point2{a[0] - b[0], a[1] - b[1]}
}

func v2Scale(a Point2, s float64) Point2 {
	return Point2{a[0] * s, a[1] * s}
}

func v2Negate(a Point2) Point2 {
	return Point2{-a[0], -a[1]}
}

func v2Dot(a, b Point2) float64 {
	return a[0]*b[0] + a[1]*b[1]
}

func v2SquaredLength(a Point2) float64 {
	return v2Dot(a, a)
}

func v2Normalize(a Point2) Point2 {
	norm := math.Sqrt(v2SquaredLength(a))
	return Point2{
		a[0] / norm,
		a[1] / norm,
	}
}
