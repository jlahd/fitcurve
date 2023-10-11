package fitcurve

import (
	"math"
	"testing"
)

func TestCurve(t *testing.T) {
	d := []Point2{
		{0.0, 0.0},
		{0.0, 0.5},
		{1.1, 1.4},
		{2.1, 1.6},
		{3.2, 1.1},
		{4.0, 0.2},
		{4.0, 0.0},
	}
	error := 4.0
	curves := FitCurve(d, error)
	if len(curves) != 1 {
		t.Errorf("Expected 1 curve, got %d", len(curves))
	}
	if curves[0][0] != d[0] {
		t.Errorf("Expected first point to be %v, got %v", d[0], curves[0][0])
	}
	if curves[0][3] != d[len(d)-1] {
		t.Errorf("Expected last point to be %v, got %v", d[len(d)-1], curves[0][len(curves[0])-1])
	}
	if curves[0][1][0] != 0 || math.Abs(curves[0][1][1]-2.1605) > 0.0001 {
		t.Errorf("Expected second point to be (0, 2.1605), got %v", curves[0][1])
	}
	if curves[0][2][0] != 4 || math.Abs(curves[0][2][1]-1.9680) > 0.0001 {
		t.Errorf("Expected third point to be (4, 1.9680), got %v", curves[0][2])
	}
}
