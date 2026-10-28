"""Independent CAD smoke check; requires cadquery-ocp (not a package dependency)."""
import math
from pathlib import Path
import sys

from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.IGESControl import IGESControl_Reader
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS


def near(actual, expected, tolerance=1e-7):
    assert abs(actual - expected) <= tolerance, (actual, expected)


def load(directory, name):
    path = directory / (name + ".igs")
    reader = IGESControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone, path
    assert reader.TransferRoots() == 1, path
    shape = reader.OneShape()
    assert BRepCheck_Analyzer(shape).IsValid(), path
    edges = TopExp_Explorer(shape, TopAbs_EDGE)
    assert edges.More(), path
    curves = []
    while edges.More():
        curves.append(BRepAdaptor_Curve(TopoDS.Edge_s(edges.Current())))
        edges.Next()
    return reader, shape, curves


def point(curve, fraction):
    u = curve.FirstParameter() + fraction * (curve.LastParameter() - curve.FirstParameter())
    p = curve.Value(u)
    return p.X(), p.Y(), p.Z()


def check_circle(curve, cx, cy, radius):
    for i in range(65):
        x, y, z = point(curve, i / 64)
        near(math.hypot(x - cx, y - cy), radius)
        near(z, 0)


def check_length(shape, expected):
    props = GProp_GProps()
    BRepGProp.LinearProperties_s(shape, props)
    near(props.Mass(), expected, 1e-5)


def main(directory):
    for name in ("line", "long_metadata"):
        reader, shape, curves = load(directory, name)
        curve, = curves
        check_length(shape, 5)
        for fraction, expected in ((0, (1, 2, 0)), (1, (4, 6, 0))):
            for actual, target in zip(point(curve, fraction), expected):
                near(actual, target)
        if name == "long_metadata":
            header = reader.IGESModel().GlobalSection()
            assert header.SendName().ToCString() == "literal,semicolon; " * 12
            assert header.AuthorName().ToCString() == "  Iges2D.jl CAD import check  "
        print(name + ": valid edge, endpoints and length checked")

    for name in ("quarter_arc", "clockwise_arc", "rational_quarter"):
        reader, shape, curves = load(directory, name)
        curve, = curves
        check_circle(curve, 0, 0, 10)
        check_length(shape, 5 * math.pi)
        for fraction, expected in ((0, (10, 0, 0)), (1, (0, 10, 0))):
            for actual, target in zip(point(curve, fraction), expected):
                near(actual, target)
        if name == "rational_quarter":
            entity = reader.IGESModel().Entity(1)
            assert entity.IsPlanar() and not entity.IsPolynomial(True)
            near(entity.UMin(), 2)
            near(entity.UMax(), 5)
        print(name + ": valid quarter circle, radius, endpoints and length checked")

    for name in ("circle", "nurbs_circle"):
        reader, shape, curves = load(directory, name)
        for curve in curves:
            check_circle(curve, 3, -2, 5)
        check_length(shape, 10 * math.pi)
        # OCCT can split type 126 at repeated knots; require a closed contour,
        # not a particular edge count. Every endpoint must have one partner.
        ends = [point(c, t) for c in curves for t in (0, 1)]
        for end in ends:
            assert sum(math.dist(end, other) < 1e-7 for other in ends) == 2
        if name == "nurbs_circle":
            assert reader.IGESModel().Entity(1).IsClosed()
        print(name + ": valid closed circle, radius and circumference checked")

    reader, shape, curves = load(directory, "polynomial")
    curve, = curves
    entity = reader.IGESModel().Entity(1)
    assert entity.IsPolynomial(True)
    near(entity.UMin(), 2)
    near(entity.UMax(), 5)
    for i in range(33):
        t = i / 32
        x, y, z = point(curve, t)
        near(x, 10 * t)
        near(y, 20 * t * (1 - t))
        near(z, 0)
    print("polynomial: valid quadratic, Bernstein values and parameter range checked")
    print("OpenCASCADE: all 8 IGES imports passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python test/cad_import.py OUTPUT_DIRECTORY")
    main(Path(sys.argv[1]))
