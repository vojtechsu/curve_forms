from utils import GF
from abc import ABC, abstractmethod


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


class NoPoint(Exception):
    pass
