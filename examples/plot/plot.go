// SPDX-License-Identifier: MPL-2.0

// Package plot is a tiny dependency-free (r, z) plotter for the Traceon examples: it maps world
// coordinates to an RGBA image and draws the things an electron-optics figure needs — electrode
// outlines, electron trajectories, a potential heatmap, and the optical axis. It exists so the
// examples can mirror upstream Traceon's matplotlib figures without pulling a plotting library
// into the add-in module.
package plot

import (
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
)

// Canvas is an (r, z) plotting surface: z runs left→right (the optical axis), r runs bottom→top.
type Canvas struct {
	img                    *image.RGBA
	w, h                   int
	zMin, zMax, rMin, rMax float64
	pad                    int
}

// New creates a Canvas of the given pixel size spanning the world box [zMin,zMax]×[rMin,rMax],
// filled with a dark background.
func New(w, h int, zMin, zMax, rMin, rMax float64) *Canvas {
	c := &Canvas{img: image.NewRGBA(image.Rect(0, 0, w, h)), w: w, h: h, zMin: zMin, zMax: zMax, rMin: rMin, rMax: rMax, pad: 28}
	bg := color.RGBA{0x1f, 0x24, 0x30, 0xff}
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c.img.SetRGBA(x, y, bg)
		}
	}
	return c
}

// px maps a world (z, r) point to a pixel coordinate (z → x, r → y, y flipped).
func (c *Canvas) px(z, r float64) (int, int) {
	inner := func(v, lo, hi, size float64) float64 {
		if hi == lo {
			return 0
		}
		return (v - lo) / (hi - lo) * size
	}
	x := float64(c.pad) + inner(z, c.zMin, c.zMax, float64(c.w-2*c.pad))
	y := float64(c.pad) + inner(r, c.rMin, c.rMax, float64(c.h-2*c.pad))
	return int(math.Round(x)), c.h - 1 - int(math.Round(y))
}

// Heatmap fills the plot area with a blue→white→red potential map: value[iz][ir] sampled on a
// uniform grid over the world box, symmetric about zero so ground reads white.
func (c *Canvas) Heatmap(values [][]float64) {
	peak := 0.0
	for _, row := range values {
		for _, v := range row {
			if a := math.Abs(v); a > peak && !math.IsInf(a, 0) && !math.IsNaN(a) {
				peak = a
			}
		}
	}
	if peak == 0 {
		peak = 1
	}
	nz, nr := len(values), 0
	if nz > 0 {
		nr = len(values[0])
	}
	for y := c.pad; y < c.h-c.pad; y++ {
		for x := c.pad; x < c.w-c.pad; x++ {
			fz := float64(x-c.pad) / float64(c.w-2*c.pad-1)
			fr := float64((c.h-1-y)-c.pad) / float64(c.h-2*c.pad-1)
			iz := clampIdx(int(fz*float64(nz-1)+0.5), nz)
			ir := clampIdx(int(fr*float64(nr-1)+0.5), nr)
			if nz == 0 || nr == 0 {
				continue
			}
			c.img.SetRGBA(x, y, divergingColor(values[iz][ir]/peak))
		}
	}
}

// Polyline draws a connected line strip through the world points (zs[i], rs[i]).
func (c *Canvas) Polyline(zs, rs []float64, col color.RGBA, width int) {
	for i := 0; i+1 < len(zs); i++ {
		x0, y0 := c.px(zs[i], rs[i])
		x1, y1 := c.px(zs[i+1], rs[i+1])
		c.line(x0, y0, x1, y1, col, width)
	}
}

// Axis draws the optical (z) axis at r = 0 as a faint horizontal line.
func (c *Canvas) Axis() {
	x0, y0 := c.px(c.zMin, 0)
	x1, y1 := c.px(c.zMax, 0)
	c.line(x0, y0, x1, y1, color.RGBA{0x3b, 0x82, 0xf6, 0x80}, 1)
}

// Save writes the canvas to a PNG file.
func (c *Canvas) Save(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := png.Encode(f, c.img); err != nil {
		_ = f.Close()
		return err
	}
	return f.Close()
}

// line draws a width-thick line between two pixels (Bresenham, with a square brush).
func (c *Canvas) line(x0, y0, x1, y1 int, col color.RGBA, width int) {
	dx, dy := abs(x1-x0), -abs(y1-y0)
	sx, sy := sign(x1-x0), sign(y1-y0)
	err := dx + dy
	for {
		c.brush(x0, y0, col, width)
		if x0 == x1 && y0 == y1 {
			return
		}
		e2 := 2 * err
		if e2 >= dy {
			err += dy
			x0 += sx
		}
		if e2 <= dx {
			err += dx
			y0 += sy
		}
	}
}

// brush paints a width×width block centred at (x, y), clipped to the image.
func (c *Canvas) brush(x, y int, col color.RGBA, width int) {
	r := width / 2
	for j := -r; j <= r; j++ {
		for i := -r; i <= r; i++ {
			px, py := x+i, y+j
			if px >= 0 && px < c.w && py >= 0 && py < c.h {
				c.img.SetRGBA(px, py, col)
			}
		}
	}
}

// divergingColor maps t∈[-1,1] to a blue→white→red colour.
func divergingColor(t float64) color.RGBA {
	if t < -1 {
		t = -1
	}
	if t > 1 {
		t = 1
	}
	lerp := func(a, b, f float64) uint8 { return uint8(a + (b-a)*f) }
	if t < 0 {
		f := t + 1 // -1→0 (blue), 0→1 (white)
		return color.RGBA{lerp(26, 242, f), lerp(51, 242, f), lerp(230, 242, f), 0xff}
	}
	return color.RGBA{lerp(242, 230, t), lerp(242, 38, t), lerp(242, 25, t), 0xff} // white→red
}

func clampIdx(i, n int) int {
	if i < 0 {
		return 0
	}
	if i >= n {
		return n - 1
	}
	return i
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func sign(a int) int {
	switch {
	case a > 0:
		return 1
	case a < 0:
		return -1
	default:
		return 0
	}
}

// Named colours used by the examples.
var (
	Electrode = color.RGBA{0xff, 0xd8, 0x33, 0xff} // amber
	Ray       = color.RGBA{0x3d, 0xdc, 0x97, 0xff} // green
	Mirror    = color.RGBA{0xc8, 0x8a, 0x4a, 0xff} // brown
)
