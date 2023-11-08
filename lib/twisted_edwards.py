from utils import GF, sqrt
from point import Point, Curve


class TwistedEdwards(Curve):
    def __init__(self, a: GF, d: GF):
        self.a = a
        self.d = d
        self.field = a.parent()
        self._form = "twiedw"
        self.params = (a, d)

    def addition(self, point1, point2):
        spoint = twiedw_sum(
            (self.a, self.d), (point1.x, point1.y), (point2.x, point2.y)
        )
        return Point(self, spoint[0], spoint[1])

    def is_infinity(self, x, y, z):
        return x == 0 and y == 1

    def infinity(self):
        return Point(self, self.field(0), self.field(1))

    def form(self):
        return "Twisted Edwards"

    def negative(self, point):
        return Point(self, -point.x, point.y, point.z)

    def check_point(self, x, y, z=1):
        return self.a * x**2 + y**2 == 1 + self.d * x**2 * y**2

    def __repr__(self):
        return f"Twisted Edwards curve {self.a}x^2+y^2=1+{self.d}x^2y^2 over F_{self.field.order()}"

    def lift_y(self, y):
        x = PolynomialRing(self.field, "x").gen()
        roots = (self.a * x**2 + y**2 - 1 - self.d * x**2 * y**2).roots()
        if roots == []:
            raise NoPoint("No such point")
        return Point(self, roots[0][0], y)

    def __iter__(self):
        try:
            for y in self.field:
                try:
                    point = self.lift_y(y)
                    yield point
                    if not point.x == 0:
                        yield -point
                except NoPoint:
                    continue
        except GeneratorExit:
            pass


def twiedw_sum(twiedw: tuple, point1: tuple, point2: tuple):
    a, d = twiedw
    x1, y1 = point1
    x2, y2 = point2
    x3 = (x1 * y2 + y1 * x2) / (1 + d * x1 * x2 * y1 * y2)
    y3 = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
    return x3, y3


def twiedw_dbl(twiedw: tuple, point: tuple):
    return twiedw_sum(twiedw, point, point)
