This is a Golang implementation of Philip J. Schneider's "Algorithm
for Automatically Fitting Digitized Curves" from the book "Graphics
Gems".

The original C code can be found at
[https://graphicsgems.org/](https://graphicsgems.org/).

## Usage

```bash
go get github.com/jlahd/fitcurve
```

```go
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
```
