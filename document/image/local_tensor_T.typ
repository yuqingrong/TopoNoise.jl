// Local tensor and blocked tensor U, drawn entirely with native Typst shapes.
// Compile with: typst compile local_tensor_T.typ local_tensor_T.pdf
#import "@preview/cetz:0.4.0": canvas, draw, tree, vector, decorations, coordinate
#import "@preview/quill:0.7.2": *
#set page(width: 760pt, height: 375pt, margin: 0pt)

#let black = rgb("111111")
#let blue = rgb("7ebede")
#let at(x, y, body) = place(top + left, dx: x, dy: y, body)
#let segment(x1, y1, x2, y2, width: 2pt) = {
  at(x1, y1, line(start: (0pt, 0pt), end: (x2 - x1, y2 - y1), stroke: (paint: black, thickness: width)))
}

#let arrowline(a, b, ratio: 20%, revert: false) = {
  draw.line(a, b, name: "line")
  if revert {
    draw.line(b, (anchor: ratio, name: "line"), mark: (end: "straight",scale: 0.55))
  } else {
    draw.line(a, (anchor: 100% - ratio, name: "line"), mark: (end: "straight",scale: 0.55))
  }
}
// Built-in Typst lines have no arrowhead option, so these native polygons
// form the small triangular heads. Their tip is at (x, y).
#let up-arrow(x, y) = at(x - 4pt, y, polygon((4pt, 0pt), (0pt, 9pt), (8pt, 9pt), fill: black))
#let right-arrow(x, y) = at(x - 9pt, y - 4pt, polygon((9pt, 4pt), (0pt, 0pt), (0pt, 8pt), fill: black))
#let northeast-arrow(x, y) = at(x - 8pt, y, polygon((8pt, 0pt), (0pt, 3pt), (3pt, 8pt), fill: black))



// Left: a local tensor T and its four directed legs.
#at(136pt, 177pt, canvas(length: 1pt, {
  arrowline((0pt, 0pt), (159pt, 0pt), ratio: 20%, revert: false)
}))

#at(109pt, 216pt, text(size: 34pt)[$T$])

// Right: a six-by-six blocked tensor U.
#for x in (490pt, 526pt, 563pt, 600pt, 640pt, 676pt) {
  segment(x, 119pt, x, 52pt)
  up-arrow(x, 51pt)
  segment(x, 254pt, x, 187pt)
  up-arrow(x, 186pt)
}
#at(470pt, 119pt, rect(width: 214pt, height: 69pt, fill: blue, stroke: (paint: rgb("003b55"), thickness: 2pt)))
#at(568pt, 135pt, text(size: 36pt)[$U$])

// Group the four physical and two virtual legs.
#at(486pt, 258pt, text(size: 26pt)[$underbrace(quad quad)_("phy")$])
#at(637pt, 258pt, text(size: 26pt)[$underbrace(quad)_("virt")$])
