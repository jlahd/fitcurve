package main

import (
	"fmt"

	"github.com/jlahd/fitcurve"
)

func main() {
	d := []fitcurve.Point2{
		{0.0, 0.0},
		{0.0, 0.5},
		{1.1, 1.4},
		{2.1, 1.6},
		{3.2, 1.1},
		{4.0, 0.2},
		{4.0, 0.0},
	}
	error := 4.0
	curves := fitcurve.FitCurve(d, error)
	if len(curves) == 1 {
		fmt.Printf("%d curve:\n", len(curves))
	} else {
		fmt.Printf("%d curves:\n", len(curves))
	}
	for _, c := range curves {
		for _, p := range c {
			fmt.Printf(" (%.4f, %.4f)", p[0], p[1])
		}
		fmt.Println()
	}
}
