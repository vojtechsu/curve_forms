from utils import GF, NoPoint, Curve


class Point:
    def __init__(self, curve: Curve, x: GF, y: GF, z: GF = 1):
        self.x = x
        self.y = y
        self.z = z
        self.field = x.parent()
        self.curve = curve
        if not self.curve.check_point(x, y, z):
            raise NoPoint("Point is not on the curve")

    def is_infinity(self):
        return self.curve.is_infinity(self.x, self.y, self.z)

    def __add__(self, other):
        return self.curve.addition(self, other)

    def __rmul__(self, scalar):
        if scalar < 0:
            return (-scalar) * (-self)
        accumulator = self.curve.infinity()
        temp = self
        while scalar > 0:
            if scalar % 2 == 1:
                accumulator = accumulator + temp
            temp = temp + temp
            scalar >>= 1
        return accumulator

    def __eq__(self, other):
        return (
            self.curve == other.curve
            and self.x == other.x
            and self.y == other.y
            and self.z == other.z
        )

    def __repr__(self):
        return f"({self.x},{self.y},{self.z})"

    def __neg__(self):
        return self.curve.negative(self)

    def __sub__(self, other):
        return self + (-other)

    def affine(self):
        assert self.z != 0
        return Point(self.curve, self.x / self.z, self.y / self.z)
