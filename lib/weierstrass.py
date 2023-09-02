from utils import GF, Curve
from point import Point


class Weierstrass(Curve):
    def __init__(self, a: GF, b: GF):
        self.a = a
        self.b = b
        self.field = a.parent()
        self._form = "shortw"
        self.params = (a, b)

    def addition(self, point1, point2):
        if point1.is_infinity():
            return point2
        if point2.is_infinity():
            return point1
        try:
            if point1 == point2:
                spoint = shortw_dbl((self.a, self.b), (point1.x, point1.y))
            else:
                spoint = shortw_sum(
                    (self.a, self.b), (point1.x, point1.y), (point2.x, point2.y)
                )
        except ZeroDivisionError:
            return self.infinity()
        return Point(self, spoint[0], spoint[1])

    def is_infinity(self, x, y, z):
        return z == 0

    def infinity(self):
        return Point(self, self.field(0), self.field(1), self.field(0))

    def form(self):
        return "Weierstrass"

    def __repr__(self):
        return (
            f"Weierstrass curve y^2=x^3+{self.a}x+{self.b} over F_{self.field.order()}"
        )

    def negative(self, point):
        return Point(self, point.x, -point.y, point.z)

    def check_point(self, x, y, z=1):
        return y**2 * z == x**3 + self.a * x * z**2 + self.b * z**3


def shortw_sum(shortw: tuple, point1: tuple, point2: tuple):
    a, b = shortw
    x1, y1 = point1
    x2, y2 = point2
    s = (y1 - y2) / (x1 - x2)
    x3 = s**2 - x1 - x2
    y3 = -y1 + s * (x1 - x3)
    return x3, y3


def shortw_dbl(shortw: tuple, point: tuple):
    a, b = shortw
    x, y = point
    s = (3 * x**2 + a) / (2 * y)
    x3 = s**2 - 2 * x
    y3 = -y + s * (x - x3)
    return x3, y3
