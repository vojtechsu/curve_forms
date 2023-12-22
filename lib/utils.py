import sage.all as sage

GF = sage.GF
PolynomialRing = sage.PolynomialRing

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
    raise Exception("The curve does not support this form.")


def sqrt(x):
    p = x.parent().order()
    sqrt_x = sage.sqrt(x)
    if not sqrt_x in x.parent():
        raise Exception("No squareroot")
    return sqrt_x if int(sqrt_x) < int(p - sqrt_x) else -sqrt_x
