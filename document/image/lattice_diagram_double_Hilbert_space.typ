// Lattice tensor-network motif, drawn entirely with native Typst elements.
// Compile with: typst compile lattice_diagram.typ lattice_diagram.pdf
#set page(width: 200pt, height: 200pt, margin: 0pt)

#let black = rgb("111111")
#let red = rgb("f3131b")
#let at(x, y, body) = place(top + left, dx: x, dy: y, body)
#let dot(x, y) = at(x - 2pt, y - 2pt, circle(radius: 2pt, fill: red))

#let segment(x1, y1, x2, y2, color: black, width: 3pt) = {
  at(x1, y1, line(start: (0pt, 0pt), end: (x2 - x1, y2 - y1), stroke: (paint: color, thickness: width)))
}

#box(width: 200pt, height: 200pt, clip: true)[
  // Black square-lattice links.
  #for x in (50pt, 100pt, 150pt) { segment(x, 10pt, x, 190pt) }
  #for y in (50pt, 100pt, 150pt) { segment(10pt, y, 190pt, y) }

  // A four-leg red tensor at each lattice site.
  #for x in (50pt, 100pt, 150pt) {
    for y in (50pt, 100pt, 150pt) {
     segment(x , y - 12.5pt, x+7.5pt, y -12.5pt-7.5pt, color: red, width: 1.5pt)
     segment(x , y + 12.5pt, x+7.5pt, y +12.5pt-7.5pt, color: red, width: 1.5pt)
    segment(x -12.5pt, y, x -12.5pt+7.5pt, y -7.5pt, color: red, width: 1.5pt)
    segment(x +12.5pt, y, x +12.5pt+7.5pt, y -7.5pt, color: red, width: 1.5pt)

     dot(x, y+12.5pt)
     dot(x, y -12.5pt)
     dot(x -12.5pt, y)
     dot(x +12.5pt, y)
    }
  }

  // Tensor label at the upper-right site.
  #at(155pt, 15pt, text(size: 18pt)[$T^(i j k ell)_(alpha beta gamma delta)$])
]
