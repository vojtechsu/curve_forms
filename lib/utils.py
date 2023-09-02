import sage.all as sage
from abc import ABC, abstractmethod


GF = sage.GF


def shortw_alpha_s_finder(shortw: tuple, all=False):
    a, b = shortw
    result = []
    z = sage.PolynomialRing(a.parent(), "z").gen()
    for alpha, _ in (z**3 + a * z + b).roots():
        if (3 * alpha**2 + a).is_square():
            s = 1 / (sqrt(3 * alpha**2 + a))
            if not all:
                return alpha, s
            result.append((alpha, s))
    if result:
        return result
    raise Exception("No Montgomery form for this curve")


def sqrt(x):
    p = x.parent().order()
    sqrt_x = sage.sqrt(x)
    if not sqrt_x in x.parent():
        raise Exception("No squareroot")
    return sqrt_x if int(sqrt_x) < int(p - sqrt_x) else -sqrt_x


class Curve:

    params = None
    field = None

    @abstractmethod
    def addition(self, point1, point2):
        pass

    @abstractmethod
    def is_infinity(self, x, y, z):
        pass

    @abstractmethod
    def infinity(self):
        pass

    @abstractmethod
    def form(self):
        pass

    def point(self, x, y):
        return Point(self, x, y)

    @abstractmethod
    def negative(self, x, y, z):
        pass

    @abstractmethod
    def check_point(self, x, y, z=1):
        pass

    @abstractmethod
    def __iter__(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        return self.field == other.field and self.params == other.params


class NoPoint(Exception):
    pass
