from utils import GF, Curve
from point import Point


class Edwards(Curve):
    def __init__(self, c: GF, d: GF):
        self.c = c
        self.d = d
        self.field = c.parent()
        self._form = "edward"
        self.params = (c, d)

    def addition(self, point1, point2):
        spoint = edward_sum(
            (self.c, self.d), (point1.x, point1.y), (point2.x, point2.y)
        )
        return Point(self, spoint[0], spoint[1])

    def is_infinity(self, x, y, z):
        return x == 0 and y == self.c

    def infinity(self):
        return Point(self, self.field(0), self.field(self.c))

    def form(self):
        return "Edwards"

    def negative(self, point):
        return Point(self, -point.x, point.y, point.z)

    def __repr__(self):
        return f"Edwards curve x^2+y^2={self.c}^2(1+{self.d}x^2y^2) over F_{self.field.order()}"

    def check_point(self, x, y, z=1):
        return x**2 + y**2 == self.c**2 * (1 + self.d * x**2 * y**2)

    def lift_y(self, y):
        x = PolynomialRing(self.field, "x").gen()
        roots = (x**2 + y**2 - self.c**2 * (1 + self.d * x**2 * y**2)).roots()
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


def edward_sum(edward: tuple, point1: tuple, point2: tuple):
    c, d = edward
    x1, y1 = point1
    x2, y2 = point2
    x3 = (x1 * y2 + y1 * x2) / (c * (1 + d * x1 * x2 * y1 * y2))
    y3 = (y1 * y2 - x1 * x2) / (c * (1 - d * x1 * x2 * y1 * y2))
    return x3, y3


def edward_dbl(edward: tuple, point: tuple):
    return edward_sum(edward, point, point)
