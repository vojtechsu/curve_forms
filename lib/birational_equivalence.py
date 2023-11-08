from utils import sqrt, shortw_alpha_s_finder
from edwards import Edwards
from weierstrass import Weierstrass
from twisted_edwards import TwistedEdwards
from montgomery import Montgomery
import point as pt


class BirationalEquivalence:
    def __init__(self, domain: pt.Curve, codomain: pt.Curve):
        self.domain = domain
        self.codomain = codomain
        self.domain_form = domain._form
        self.codomain_form = codomain._form
        self.field = self.domain.field

        module = __import__(__name__)
        self._mapping_function = getattr(
            module, f"{self.domain_form}_to_{self.codomain_form}_point"
        )

    @classmethod
    def to_edwards(cls, domain: pt.Curve):
        source_form = domain._form
        assert source_form != "edward"
        module = __import__(__name__)
        bir_eq = getattr(module, f"{source_form}_to_edward")
        edward = Edwards(*bir_eq(domain.params))
        return cls(domain, edward)

    @classmethod
    def to_weierstrass(cls, domain: pt.Curve):
        source_form = domain._form
        assert source_form != "shortw"
        module = __import__(__name__)
        bir_eq = getattr(module, f"{source_form}_to_shortw")
        shortw = Weierstrass(*bir_eq(domain.params))
        return cls(domain, shortw)

    @classmethod
    def to_twisted_edwards(cls, domain: pt.Curve):
        source_form = domain._form
        assert source_form != "twiedw"
        module = __import__(__name__)
        bir_eq = getattr(module, f"{source_form}_to_twiedw")
        twiedw = TwistedEdwards(*bir_eq(domain.params))
        return cls(domain, twiedw)

    @classmethod
    def to_montgomery(cls, domain: pt.Curve):
        source_form = domain._form
        assert source_form != "montgo"
        module = __import__(__name__)
        bir_eq = getattr(module, f"{source_form}_to_montgo")
        montgo = Montgomery(*bir_eq(domain.params))
        return cls(domain, montgo)

    def __call__(self, point: pt.Point):
        point = point.affine()
        x, y = point.x, point.y
        x, y = self._mapping_function(self.domain.params, (x, y), self.codomain.params)
        return pt.Point(self.codomain, x, y, self.field(1))


def shortw_to_montgo(shortw: tuple):
    alpha, s = shortw_alpha_s_finder(shortw)
    montgo_a, montgo_b = 3 * alpha * s, s
    return montgo_a, montgo_b


def shortw_to_montgo_point(shortw: tuple, point: tuple, montgo: tuple):
    x, y = point
    a, b = montgo
    alpha, s = a / (3 * b), b
    return s * (x - alpha), s * y


def montgo_to_shortw(montgo: tuple):
    a, b = montgo
    shortw_a, shortw_b = (3 - a**2) / (3 * b**2), (2 * a**3 - 9 * a) / (
        27 * b**3
    )
    return shortw_a, shortw_b


def montgo_to_shortw_point(montgo: tuple, point: tuple, shortw: tuple):
    x, y = point
    a, b = montgo
    return (3 * x + a) / (3 * b), y / b


def shortw_to_twiedw(shortw: tuple):
    a, b = shortw
    alpha, s = shortw_alpha_s_finder(shortw)
    twiedw_a, twiedw_d = 3 * alpha + 2 / s, 3 * alpha - 2 / s
    return twiedw_a, twiedw_d


def shortw_to_twiedw_point(shortw: tuple, point: tuple, twiedw: tuple):
    x, y = point
    a, d = twiedw
    alpha, s = (a + d) / 6, 4 / (a - d)
    return (x - alpha) / y, (s * (x - alpha) - 1) / (s * (x - alpha) + 1)


def twiedw_to_shortw(twiedw):
    a, d = twiedw
    shortw_a = -(a**2 + 14 * d * a + d**2) / 48
    shortw_b = (a + d) * (-(a**2) + 34 * a * d - d**2) / 864
    return shortw_a, shortw_b


def twiedw_to_shortw_point(twiedw: tuple, point: tuple, shortw: tuple):
    x, y = point
    a, d = twiedw
    u = (5 * a + a * y - 5 * d * y - d) / (12 - 12 * y)
    v = (a + a * y - d * y - d) / (4 * x - 4 * x * y)
    return u, v


def montgo_to_twiedw(montgo: tuple):
    a, b = montgo
    twiedw_a, twiedw_d = (a + 2) / b, (a - 2) / b
    return twiedw_a, twiedw_d


def montgo_to_twiedw_point(montgo: tuple, point: tuple, twiedw: tuple):
    x, y = point
    return x / y, (x - 1) / (x + 1)


def twiedw_to_montgo(twiedw: tuple):
    a, d = twiedw
    montgo_a, montgo_b = 2 * (a + d) / (a - d), 4 / (a - d)
    return montgo_a, montgo_b


def twiedw_to_montgo_point(twiedw: tuple, point: tuple, montgo: tuple):
    x, y = point
    return (1 + y) / (1 - y), (1 + y) / ((1 - y) * x)


def twiedw_to_edward(twiedw: tuple):
    a, d = twiedw
    sqrt_a = sqrt(1 / a)
    edward_a, edward_d = sqrt_a, d * a
    return edward_a, edward_d


def twiedw_to_edward_point(twiedw: tuple, point: tuple, edward: tuple):
    x, y = point
    c, d = edward
    return x, y * c


def edward_to_twiedw(edward: tuple):
    c, d = edward
    twiedw_a, twiedw_d = 1 / c**2, c**2 * d
    return twiedw_a, twiedw_d


def edward_to_twiedw_point(edward: tuple, point: tuple, twiedw: tuple):
    c, d = edward
    x, y = point
    return x, y / c


def edward_to_shortw(edward: tuple):
    c, d = edward
    shortw_a, shortw_b = -(1 + 14 * d * c**4 + c**8 * d**2) / (48 * c**4), (
        1 + c**4 * d
    ) * (-1 + 34 * d * c**4 - c**8 * d**2) / (864 * c**6)
    return shortw_a, shortw_b


def edward_to_shortw_point(edward: tuple, point: tuple, shortw: tuple):
    c, d = edward
    x, y = point
    u, v = (5 * c + y - 5 * c**4 * d * y - c**5 * d) / (
        12 * c**3 - 12 * y * c**2
    ), (c + y - d * y * c**4 - c**5 * d) / (4 * x * c**3 - 4 * x * y * c**2)
    return u, v


def shortw_to_edward(shortw: tuple):
    for alpha, s in shortw_alpha_s_finder(shortw, all=True):
        try:
            t = sqrt(s / (3 * s * alpha + 2))
        except:
            continue
        a, b = shortw
        edward_c, edward_d = t, -4 * a - 3 * alpha**2
        return edward_c, edward_d


def shortw_to_edward_point(shortw: tuple, point: tuple, edward: tuple):
    x, y = point
    a, b = shortw
    c, d = edward
    s = 1 / sqrt(-3 * a - d)
    alpha = (1 / c**2 - 2 / s) / 3
    u, v = (x - alpha) / y, (s * (x - alpha) - 1) / (s * (x - alpha) + 1) * c
    return u, v


def montgo_to_edward(montgo: tuple):
    a, b = montgo
    edward_c, edward_d = sqrt(b / (a + 2)), (a**2 - 4) / (b**2)
    return edward_c, edward_d


def montgo_to_edward_point(montgo: tuple, point: tuple, edward: tuple):
    a, b = montgo
    x, y = point
    return x / y, (x - 1) / (x + 1) * sqrt(b / (a + 2))


def edward_to_montgo(edward: tuple):
    c, d = edward
    montgo_a, montgo_b = (2 + 2 * c**4 * d) / (1 - c**4 * d), (4 * c**2) / (
        1 - c**4 * d
    )
    return montgo_a, montgo_b


def edward_to_montgo_point(edward: tuple, point: tuple, montgo: tuple):
    c, d = edward
    x, y = point
    return (c + y) / (c - y), (c + y) / ((c - y) * x)
