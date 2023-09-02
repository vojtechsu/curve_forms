from sage.all import EllipticCurve, GF, PolynomialRing
import weierstrass
import edwards
import twisted_edwards
import montgomery
from point import Point
import birationa_equivalence as be


class TestTransformations:

    F = GF(1009)
    w_vec = {
        "a": F(866),
        "b": F(208),
        "P": (F(353), F(449)),
        "Q": (F(924), F(356)),
        "R": (F(542), F(665)),
    }
    m_vec = {
        "a": F(32),
        "b": F(733),
        "P": (F(98), F(183)),
        "Q": (F(915), F(626)),
        "R": (F(402), F(98)),
    }
    tw_vec = {
        "a": F(519),
        "d": F(636),
        "P": (F(188), F(480)),
        "Q": (F(593), F(435)),
        "R": (F(766), F(677)),
    }
    ed_vec = {
        "c": F(480),
        "d": F(141),
        "P": (F(188), F(348)),
        "Q": (F(593), F(946)),
        "R": (F(766), F(62)),
    }
    E = EllipticCurve(F, [w_vec["a"], w_vec["b"]])
    P = E(*w_vec["P"])
    Q = E(*w_vec["Q"])
    R = E(*w_vec["R"])
    shortw = w_vec["a"], w_vec["b"]
    edward = ed_vec["c"], ed_vec["d"]
    twiedw = tw_vec["a"], tw_vec["d"]
    montgo = m_vec["a"], m_vec["b"]
    assert P + Q == R

    def test_edward(self):
        formula_test_help(
            self.P,
            self.Q,
            self.R,
            be.shortw_to_edward,
            be.shortw_to_edward_point,
            edwards.edward_sum,
        )

    def test_twiedw(self):
        formula_test_help(
            self.P,
            self.Q,
            self.R,
            be.shortw_to_twiedw,
            be.shortw_to_twiedw_point,
            twisted_edwards.twiedw_sum,
        )

    def test_montgo(self):
        formula_test_help(
            self.P,
            self.Q,
            self.R,
            be.shortw_to_montgo,
            be.shortw_to_montgo_point,
            montgomery.montgo_sum,
        )

    def test_input_points(self):
        x, y = PolynomialRing(self.F, ["x", "y"]).gens()
        for point in ["P", "Q", "R"]:
            x, y = self.w_vec[point]
            assert y**2 == x**3 + self.w_vec["a"] * x + self.w_vec["b"]
            x, y = self.m_vec[point]
            assert self.m_vec["b"] * y**2 == x**3 + self.m_vec["a"] * x**2 + x
            x, y = self.tw_vec[point]
            assert (
                self.tw_vec["a"] * x**2 + y**2
                == 1 + self.tw_vec["d"] * x**2 * y**2
            )
            x, y = self.ed_vec[point]
            assert x**2 + y**2 == self.ed_vec["c"] ** 2 * (
                1 + self.ed_vec["d"] * x**2 * y**2
            )
        assert weierstrass.shortw_sum(self.shortw, self.P[:2], self.Q[:2]) == self.R[:2]

    def test_shortw_to_other(self):
        assert self.montgo == be.shortw_to_montgo(self.shortw)
        assert self.twiedw == be.shortw_to_twiedw(self.shortw)
        assert self.edward == be.shortw_to_edward(self.shortw)

    def test_montgo_to_other(self):
        assert self.shortw == be.montgo_to_shortw(self.montgo)
        assert self.twiedw == be.montgo_to_twiedw(self.montgo)
        assert self.edward == be.montgo_to_edward(self.montgo)

    def test_edward_to_other(self):
        assert self.shortw == be.edward_to_shortw(self.edward)
        assert self.twiedw == be.edward_to_twiedw(self.edward)
        assert self.montgo == be.edward_to_montgo(self.edward)

    def test_twiedw_to_other(self):
        assert self.shortw == be.twiedw_to_shortw(self.twiedw)
        assert self.montgo == be.twiedw_to_montgo(self.twiedw)
        assert self.edward == be.twiedw_to_edward(self.twiedw)

    def test_shortw_to_other_point(self):
        for point in ["P", "Q", "R"]:
            wpoint = self.w_vec[point]
            assert self.m_vec[point] == be.shortw_to_montgo_point(
                self.shortw, wpoint, self.montgo
            )
            assert self.tw_vec[point] == be.shortw_to_twiedw_point(
                self.shortw, wpoint, self.twiedw
            )
            assert self.ed_vec[point] == be.shortw_to_edward_point(
                self.shortw, wpoint, self.edward
            )

    def test_montgo_to_other_point(self):
        for point in ["P", "Q", "R"]:
            mpoint = self.m_vec[point]
            assert self.w_vec[point] == be.montgo_to_shortw_point(
                self.montgo, mpoint, self.shortw
            )
            assert self.tw_vec[point] == be.montgo_to_twiedw_point(
                self.montgo, mpoint, self.twiedw
            )
            assert self.ed_vec[point] == be.montgo_to_edward_point(
                self.montgo, mpoint, self.edward
            )

    def test_edward_to_other_point(self):
        for point in ["P", "Q", "R"]:
            epoint = self.ed_vec[point]
            assert self.m_vec[point] == be.edward_to_montgo_point(
                self.edward, epoint, self.montgo
            )
            assert self.tw_vec[point] == be.edward_to_twiedw_point(
                self.edward, epoint, self.twiedw
            )
            assert self.w_vec[point] == be.edward_to_shortw_point(
                self.edward, epoint, self.shortw
            )

    def test_twiedw_to_other_point(self):
        for point in ["P", "Q", "R"]:
            tpoint = self.tw_vec[point]
            assert self.m_vec[point] == be.twiedw_to_montgo_point(
                self.twiedw, tpoint, self.montgo
            )
            assert self.w_vec[point] == be.twiedw_to_shortw_point(
                self.twiedw, tpoint, self.shortw
            )
            assert self.ed_vec[point] == be.twiedw_to_edward_point(
                self.twiedw, tpoint, self.edward
            )


def formula_test_help(P, Q, R, transformation, point_mapping, summation):
    shortw = P.curve().a4(), P.curve().a6()
    otherf = transformation(shortw)
    mP = point_mapping(shortw, P[:2], otherf)
    mQ = point_mapping(shortw, Q[:2], otherf)
    mR = point_mapping(shortw, R[:2], otherf)
    assert summation(otherf, mP, mQ) == mR, {f"{summation(shortw,mP,mQ)}, {mR}"}


def test_weierstrass():
    F = GF(101)
    a, b = F(1), F(2)
    W = weierstrass.Weierstrass(a, b)
    I = W.infinity()
    x, y, z = I.x, I.y, I.z
    assert (x, y, z) == (0, 1, 0)
    assert I.is_infinity()
    P = Point(W, F(77), F(30))
    P6 = 6 * P
    assert (P6.x, P6.y, P6.z) == (F(18), F(14), F(1))
    Q = Point(W, F(17), F(36))
    R = P + Q
    assert (R.x, R.y, R.z) == (F(6), F(74), F(1))
    assert (P + I) == P

    WM = be.BirationalEquivalence.to_montgomery(W)
    WE = be.BirationalEquivalence.to_edwards(W)
    WT = be.BirationalEquivalence.to_twisted_edwards(W)
    assert WE(P) + WE(Q) == WE(R)
    assert WM(P) + WM(Q) == WM(R)
    assert WT(P) + WT(Q) == WT(R)


def test_montgomery():
    F = GF(101)
    a, b = F(49), F(51)
    M = montgomery.Montgomery(a, b)
    I = M.infinity()
    x, y, z = I.x, I.y, I.z
    assert (x, y, z) == (0, 1, 0)
    assert I.is_infinity()
    P = Point(M, F(39), F(15))
    P6 = 6 * P
    assert (P6.x, P6.y, P6.z) == (F(60), F(7), F(1))
    Q = Point(M, F(9), F(18))
    R = P + Q
    assert (R.x, R.y, R.z) == (F(54), F(37), F(1))
    assert (P + I) == P

    tW = be.BirationalEquivalence.to_weierstrass(M)
    tE = be.BirationalEquivalence.to_edwards(M)
    tT = be.BirationalEquivalence.to_twisted_edwards(M)
    assert tE(P) + tE(Q) == tE(R)
    assert tW(P) + tW(Q) == tW(R)
    assert tT(P) + tT(Q) == tT(R)


def test_twistededwards():
    F = GF(101)
    a, d = F(1), F(94)
    T = twisted_edwards.TwistedEdwards(a, d)
    I = T.infinity()
    x, y, z = I.x, I.y, I.z
    assert (x, y) == (0, 1)
    assert I.is_infinity()
    P = Point(T, F(51), F(21))
    P6 = 6 * P
    assert (P6.x, P6.y, P6.z) == (F(39), F(13), F(1))
    Q = Point(T, F(43), F(6))
    R = P + Q
    assert (R.x, R.y, R.z) == (F(97), F(23), F(1))
    assert (P + I) == P

    tW = be.BirationalEquivalence.to_weierstrass(T)
    tE = be.BirationalEquivalence.to_edwards(T)
    tM = be.BirationalEquivalence.to_montgomery(T)
    assert tE(P) + tE(Q) == tE(R)
    assert tW(P) + tW(Q) == tW(R)
    assert tM(P) + tM(Q) == tM(R)


def test_edwards():
    F = GF(101)
    c, d = F(36), F(55)
    T = edwards.Edwards(c, d)
    I = T.infinity()
    x, y, z = I.x, I.y, I.z
    assert (x, y) == (0, c)
    assert I.is_infinity()
    P = Point(T, F(29), F(100))
    P6 = 6 * P
    assert (P6.x, P6.y, P6.z) == (F(93), F(85), F(1))
    Q = Point(T, F(46), F(17))
    R = P + Q
    assert (R.x, R.y, R.z) == (F(10), F(80), F(1))
    assert (P + I) == P

    tW = be.BirationalEquivalence.to_weierstrass(domain=T)
    tE = be.BirationalEquivalence.to_twisted_edwards(T)
    tT = be.BirationalEquivalence.to_montgomery(T)
    assert tE(P) + tE(Q) == tE(R)
    assert tW(P) + tW(Q) == tW(R)
    assert tT(P) + tT(Q) == tT(R)
