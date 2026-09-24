#import "@preview/cetz:0.5.2"

// Reusable site and gate primitives for a two-dimensional circuit.
// The canvas is tall enough for the six gate-offset levels above the lattice.
#let panel-width = 493
#let panel-height = 400

#let foreground = rgb("#2f3440")
#let circuit-stroke = (
  paint: foreground,
  thickness: 1.5pt,
  cap: "round",
  join: "round",
)
#let gate-offset = 30

// Draw a named site anchor. Gates and wires can refer to `name` later.
#let port(name, x, y, fill-color: white) = {
  cetz.draw.circle(
    (x, y),
    name: name,
    radius: 8.5,
    fill: fill-color,
    stroke: circuit-stroke,
  )
}

// Draw an H gate above a named site instead of covering the site marker.
// With this file's 1pt canvas unit, the default offset is 30pt.
#let hadamard(
  site,
  offset: gate-offset,
  size: 20,
  fill-color: white,
) = {
  cetz.draw.rect(
    (rel: (-size / 2, offset - size / 2), to: site),
    (rel: (size / 2, offset + size / 2), to: site),
    radius: 2,
    fill: fill-color,
    stroke: circuit-stroke,
  )
  cetz.draw.content(
    (rel: (0, offset), to: site),
    text(size: 15pt, fill: foreground)[$H$],
  )
}

// Draw a CNOT above two named sites. The first site is the control and the
// second is the target; both symbols receive the same upward offset.
#let cnot(
  control-site,
  target-site,
  offset: gate-offset,
  target-radius: 8.5,
  fill-color: white,
) = {
  let control = (rel: (0, offset), to: control-site)
  let target = (rel: (0, offset), to: target-site)

  cetz.draw.line(control, target, stroke: circuit-stroke)
  cetz.draw.circle(
    control,
    radius: 4,
    fill: foreground,
    stroke: none,
  )
  cetz.draw.circle(
    target,
    radius: target-radius,
    fill: fill-color,
    stroke: circuit-stroke,
  )
  cetz.draw.line(
    (rel: (-target-radius, offset), to: target-site),
    (rel: (target-radius, offset), to: target-site),
    stroke: circuit-stroke,
  )
  cetz.draw.line(
    (rel: (0, offset - target-radius), to: target-site),
    (rel: (0, offset + target-radius), to: target-site),
    stroke: circuit-stroke,
  )
}

#set page(
  width: panel-width * 1pt,
  height: panel-height * 1pt,
  margin: 0pt,
)
#set text(
  font: "New Computer Modern",
  fill: foreground,
)

#cetz.canvas(length: 1pt, {
  import cetz.draw: *

  // Establish the full panel bounds without painting a background.
  rect((0, 0), (panel-width, panel-height), fill: none, stroke: none)

  // Nearest-neighbor edges of the sheared 3 x 3 square lattice. Drawing
  // these first lets the opaque site markers hide each line's center end.
  for (from, to) in (
    ((66, 20), (177, 20)),
    ((177, 20), (287, 20)),
    ((105, 107), (215, 107)),
    ((215, 107), (326, 107)),
    ((146, 196), (256, 196)),
    ((256, 196), (366, 196)),
    ((66, 20), (105, 107)),
    ((105, 107), (146, 196)),
    ((177, 20), (215, 107)),
    ((215, 107), (256, 196)),
    ((287, 20), (326, 107)),
    ((326, 107), (366, 196)),
  ) {
    line(from, to, stroke: circuit-stroke)
  }

  // A short wire above every site provides a mounting point for a gate.
  for (x, y) in (
    (66, 20),
    (177, 20),
    (287, 20),
    (105, 107),
    (215, 107),
    (326, 107),
    (146, 196),
    (256, 196),
    (366, 196),
  ) {
    line((x, y), (x, y + 6*gate-offset), stroke: circuit-stroke)
  }

  // Site anchors, indexed as site(x)(y) from the lower-left corner.
  for (name, x, y) in (
    ("site11", 66, 20),
    ("site21", 177, 20),
    ("site31", 287, 20),
    ("site12", 105, 107),
    ("site22", 215, 107),
    ("site32", 326, 107),
    ("site13", 146, 196),
    ("site23", 256, 196),
    ("site33", 366, 196),
  ) {
    port(name, x, y)
  }


  hadamard("site21")
  cnot("site11","site21", offset: 2*gate-offset)

  hadamard("site22")
  cnot("site12","site22",offset: 2*gate-offset)
  cnot("site13","site22",offset: 3*gate-offset)
  cnot("site23","site22",offset: 4*gate-offset)

hadamard("site33", offset: 3.5*gate-offset)
cnot("site23","site33", offset: 4.5*gate-offset )

hadamard("site31", offset: 3*gate-offset)
cnot("site21","site31",offset: 4*gate-offset)
cnot("site22","site31",offset: 5*gate-offset)
cnot("site32","site31",offset: 5.8*gate-offset)
 
port("site21", 177, 20, fill-color: red)
port("site22", 215, 107, fill-color: red)
port("site31", 287, 20, fill-color: red)
port("site33", 366, 196, fill-color: red)
})