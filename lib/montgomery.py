from utils import GF, Curve
import point as pt


class Montgomery(Curve):
    def __init__(self, a: GF, b: GF):
        self.a = a
        self.b = b
        self.field = a.parent()
        self._form = "montgo"
        self.params = (a, b)

    def addition(self, point1, point2):
        if point1.is_infinity():
            return point2
        if point2.is_infinity():
            return point1
        try:
            if point1 == point2:
                spoint = montgo_dbl((self.a, self.b), (point1.x, point1.y))
            else:
                spoint = montgo_sum(
                    (self.a, self.b), (point1.x, point1.y), (point2.x, point2.y)
                )
        except ZeroDivisionError:
            return self.infinity()
        return pt.Point(self, spoint[0], spoint[1])

    def is_infinity(self, x, y, z):
        return z == 0

    def infinity(self):
        return pt.Point(self, self.field(0), self.field(1), self.field(0))

    def form(self):
        return "Montgomery"

    def negative(self, point):
        return pt.Point(self, point.x, -point.y, point.z)

    def check_point(self, x, y, z=1):
        return self.b * y**2 * z == x**3 + self.a * x**2 * z + x * z**2


def montgo_sum(montgo: tuple, point1: tuple, point2: tuple):
    a, b = montgo
    x1, y1 = point1
    x2, y2 = point2
    x3 = b * (y2 - y1) ** 2 / (x2 - x1) ** 2 - a - x1 - x2
    y3 = (
        (2 * x1 + x2 + a) * (y2 - y1) / (x2 - x1)
        - b * (y2 - y1) ** 3 / (x2 - x1) ** 3
        - y1
    )
    return x3, y3


def montgo_dbl(montgo: tuple, point: tuple):
    a, b = montgo
    x1, y1 = point
    x3 = b * (3 * x1**2 + 2 * a * x1 + 1) ** 2 / ((2 * b * y1) ** 2) - a - x1 - x1
    y3 = (
        (3 * x1 + a) * (3 * x1**2 + 2 * a * x1 + 1) / (2 * b * y1)
        - b * (3 * x1**2 + 2 * a * x1 + 1) ** 3 / ((2 * b * y1) ** 3)
        - y1
    )
    return x3, y3
