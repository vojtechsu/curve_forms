from sage.all import EllipticCurve, GF, ZZ, PolynomialRing, sqrt, EllipticCurve_from_j
from abc import ABC, abstractmethod

def weier_to_mont(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    for alpha,_ in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            A, B = 3*alpha*s, s
            u,v = s*(x-alpha),s*y
            return A,B,u,v

def mont_to_weier(a,b,x,y):
    A, B = (3-a**2)/(3*b**2), (2*a**3-9*a)/(27*b**3)
    u,v = (3*x+a)/(3*b),y/b
    return A,B,u,v

def weier_to_twedw(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    for alpha,_ in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            A = 3*alpha+2/s
            D = 3*alpha-2/s
            u,v = (x-alpha)/y,(s*(x-alpha)-1)/(s*(x-alpha)+1)
            return A,D,u,v

def twedw_to_weier(a,d,x,y):
    A = -(a**2+14*d*a+d**2)/48
    B = (a+d)*(-a**2+34*a*d-d**2)/864
    u = (5*a+a*y-5*d*y-d)/(12-12*y)
    v = (a+a*y-d*y-d)/(4*x-4*x*y)
    return A,B,u,v

def mont_to_twedw(a,b,x,y):
    A,D = (a+2)/b,(a-2)/b
    u,v = x/y,(x-1)/(x+1)
    return A,D,u,v

def twedw_to_mont(a,d,x,y):
    A,B = 2*(a+d)/(a-d),4/(a-d)
    u,v = (1+y)/(1-y),(1+y)/((1-y)*x)
    return A,B,u,v

def twedw_to_edwar(a,d,x,y):
    p = a.parent().order()
    sa = sqrt(1/a) if ZZ(sqrt(1/a))<ZZ(p-sqrt(1/a)) else -sqrt(1/a)
    C,D = sa,d*a
    u,v = x,y*C
    return C,D,u,v

def edwar_to_twedw(c,d,x,y):
    u,v = x,y/c
    A,D = 1/c**2, c**2*d
    return A,D,u,v

def edwar_to_weier(c,d,x,y):
    u,v = (5*c+y-5*c**4*d*y-c**5*d)/(12*c**3-12*y*c**2),(c+y-d*y*c**4-c**5*d)/(4*x*c**3-4*x*y*c**2)
    A,B = -(1+14*d*c**4+c**8*d**2)/(48*c**4), (1+c**4*d)*(-1+34*d*c**4-c**8*d**2)/(864*c**6)
    return A,B,u,v


def weier_to_edwar(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    p = a.parent().order()
    for alpha,ex in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            if not (s/(3*s*alpha+2)).is_square():
                continue
            t = sqrt(s/(3*s*alpha+2))
            t = t if ZZ(t)<ZZ(p-t) else p-t
            C, D = t,-4*a-3*alpha**2
            u,v = (x-alpha)/y,(s*(x-alpha)-1)/(s*(x-alpha)+1)*t
            return C,D,u,v

def mont_to_edwar(a,b,x,y):
    C,D = sqrt(b/(a+2)),(a**2-4)/(b**2)
    u,v = x/y,(x-1)/(x+1)*sqrt(b/(a+2))
    return C,D,u,v

def edwar_to_mont(c,d,x,y):
    A,B = (2+2*c**4*d)/(1-c**4*d),(4*c**2)/(1-c**4*d)
    u,v = (c+y)/(c-y),(c+y)/((c-y)*x)
    return A,B,u,v


def edwar_sum(c,d,x1,y1,x2,y2):
    x3 = (x1*y2+y1*x2)/(c*(1+d*x1*x2*y1*y2))
    y3 = (y1*y2-x1*x2)/(c*(1-d*x1*x2*y1*y2))
    return x3,y3

def edwar_dbl(c,d,x,y):
    return edw_sum(c,d,x,y,x,y)

def mont_sum(a,b,x1,y1,x2,y2):
    x3 = b*(y2-y1)**2/(x2-x1)**2-a-x1-x2
    y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)**3/(x2-x1)**3-y1
    return x3,y3

def mont_dbl(a,b,x1,y1):
    x3 = b*(3*x1**2+2*a*x1+1)**2/((2*b*y1)**2)-a-x1-x1
    y3 = (3*x1+a)*(3*x1**2+2*a*x1+1)/(2*b*y1)-b*(3*x1**2+2*a*x1+1)**3/((2*b*y1)**3)-y1
    return x3,y3

def twedw_sum(a,d,x1,y1,x2,y2):
    x3 = (x1*y2+y1*x2)/(1+d*x1*x2*y1*y2)
    y3 = (y1*y2-a*x1*x2)/(1-d*x1*x2*y1*y2)
    return x3,y3

def twedw_dbl(a,b,x,y):
    return twedw_sum(a,b,x,y,x,y)


def test_formula(P,Q,R,transformation,summation):
    a,b = P.curve().a4(), P.curve().a6()
    s,t,x1,y1 = transformation(a,b,P[0],P[1])
    _,_,x2,y2 = transformation(a,b,Q[0],Q[1])
    _,_,x3,y3 = transformation(a,b,R[0],R[1])
    assert summation(s,t,x1,y1,x2,y2)==(x3,y3)

def test_edwar(P,Q,R):
    test_formula(P,Q,R,weier_to_edwar,edwar_sum)

def test_twedw(P,Q,R):
    test_formula(P,Q,R,weier_to_twedw,twedw_sum)

def test_mont(P,Q,R):
    test_formula(P,Q,R,weier_to_mont,mont_sum)

def find_random_edwards(F:GF, cofactor = 4):
    z = PolynomialRing(F,'z').gen()
    p = F.order()
    while True:
        E = EllipticCurve_from_j(F.random_element())
        a,b = E.a4(),E.a6()
        for alpha,_ in (z**3+a*z+b).roots():
            if (3*alpha**2+a).is_square():
                s = 1/(sqrt(3*alpha**2+a))
                if not (s/(3*s*alpha+2)).is_square():
                    continue
                order = E.order()
                if order%cofactor!=0:continue
                if not (order//cofactor).is_prime():continue
                t = sqrt(s/(3*s*alpha+2))
                t = t if ZZ(t)<ZZ(p-t) else p-t
                c, d = t,-4*a-3*alpha**2
                return Edwards(c,d)

def find_random_weierstrass(F:GF,cofactor=1,a=None,b=None):
    while True:
        if a is not None:
            try:
                E = EllipticCurve(F,[a,F.random_element()])
            except ArithmeticError:
                continue
        if b is not None:
            try:
                E = EllipticCurve(F,[F.random_element(),b])
            except ArithmeticError:
                continue
        if a is None and b is None:
            E = EllipticCurve_from_j(F.random_element())
        if cofactor is None:
            return Weierstrass(E.a4(),E.a6())
        order = E.order()
        if order%cofactor!=0:continue
        if not (order//cofactor).is_prime():continue
        return Weierstrass(E.a4(),E.a6())


def find_random_twistededwards(F: GF, cofactor = 4):
    z = PolynomialRing(F,'z').gen()
    while True:
        E = EllipticCurve_from_j(F.random_element())
        a,b = E.a4(),E.a6()
        for alpha,_ in (z**3+a*z+b).roots():
            if (3*alpha**2+a).is_square():
                s = 1/(sqrt(3*alpha**2+a))
                order = E.order()
                if order%cofactor!=0:continue
                if not (order//cofactor).is_prime():continue
                a,d = 3*alpha+2/s, 3*alpha-2/s
                return TwistedEdwards(a,d)


class Curve:
    @abstractmethod
    def addition(self,point1,point2):
        pass

    @abstractmethod
    def is_infinity(self,x,y,z):
        pass

    @abstractmethod
    def infinity(self):
        pass

    @abstractmethod
    def form(self):
        pass

    def point(self,x,y):
        return Point(self,x,y)

    @abstractmethod
    def negative(self,x,y,z):
        pass

    @abstractmethod
    def check_point(self,x,y,z=1):
        pass

    @abstractmethod
    def __iter__(self):
        pass


class NoPoint(Exception):
    pass

class Weierstrass(Curve):
    def __init__(self,a: GF,b: GF):
        self.a = a
        self.b = b
        self.field = a.parent()
        self.sage_curve = EllipticCurve(self.field,[a,b])

    def addition(self,point1,point2):
        e = self.sage_curve
        if point1.is_infinity():
            return point2
        if point2.is_infinity():
            return point1
        spoint = e(point1.x,point1.y)+e(point2.x,point2.y)
        return Point(self,spoint[0],spoint[1],spoint[2])

    def is_infinity(self,x,y,z):
        return z==0

    def infinity(self):
        return Point(self,self.field(0),self.field(1),self.field(0))

    def form(self):
        return "Weierstrass"

    def __repr__(self):
        return f"Weierstrass curve y^2=x^3+{self.a}x+{self.b} over F_{self.field.order()}"

    def negative(self,point):
        return Point(self,point.x,-point.y,point.z)

    def to_edwards(self,point=None):
        x,y = (self.field(1),self.field(1)) if point is None else (point.x,point.y)
        ed = weier_to_edwar(self.a,self.b,x,y)
        if ed is None:
            return
        c,d,x,y = ed
        edward = Edwards(c,d)
        if point is None:
            return edward
        return Point(edward,x,y)

    def check_point(self,x,y,z=1):
        return y**2*z==x**3+self.a*x*z**2+self.b*z**3

    def __iter__(self):
        for point in self.sage_curve:
            yield Point(self,point[0],point[1],point[2])



class Montgomery(Curve):
    def __init__(self,a: GF,b: GF):
        self.a = a
        self.b = b
        self.field = a.parent()

    def addition(self,point1,point2):
        if point1.is_infinity():
            return point2
        if point2.is_infinity():
            return point1
        try:
            if point1==point2:
                spoint = mont_dbl(self.a,self.b,point1.x,point1.y)
            else:
                spoint = mont_sum(self.a,self.b,point1.x,point1.y,point2.x,point2.y)
        except ZeroDivisionError:
            return self.infinity()
        return Point(self,spoint[0],spoint[1])

    def is_infinity(self,x,y,z):
        return z==0

    def infinity(self):
        return Point(self,self.field(0),self.field(1),self.field(0))

    def form(self):
        return "Montgomery"

    def negative(self,point):
        return Point(self,point.x,-point.y,point.z)

    def check_point(self,x,y,z=1):
        return self.b*y**2*z==x**3+self.a*x**2*z+x*z**2


class TwistedEdwards(Curve):
    def __init__(self,a: GF,d: GF):
        self.a = a
        self.d = d
        self.field = a.parent()

    def addition(self,point1,point2):
        spoint = twedw_sum(self.a,self.d,point1.x,point1.y,point2.x,point2.y)
        return Point(self,spoint[0],spoint[1])

    def is_infinity(self,x,y,z):
        return x==0 and y==1

    def infinity(self):
        return Point(self,self.field(0),self.field(1))

    def form(self):
        return "Twisted Edwards"

    def negative(self,point):
        return Point(self,-point.x,point.y,point.z)

    def check_point(self,x,y,z=1):
        return self.a*x**2+y**2==1+self.d*x**2*y**2

    def __repr__(self):
        return f"Twisted Edwards curve {self.a}x^2+y^2=1+{self.d}x^2y^2 over F_{self.field.order()}"


    def lift_y(self,y):
        x = PolynomialRing(self.field,'x').gen()
        roots = (self.a*x**2+y**2-1-self.d*x**2*y**2).roots()
        if roots==[]:
            raise NoPoint("No such point")
        return Point(self,roots[0][0],y)

    def __iter__(self):
        try:
            for y in self.field:
                try:
                    point = self.lift_y(y)
                    yield point
                    if not point.x==0: yield -point
                except NoPoint:
                    continue
        except GeneratorExit:
                pass

    def to_weierstrass(self,point=None):
        x,y = (self.field(1),self.field(1)) if point is None else (point.x,point.y)
        we = twedw_to_weier(self.a,self.d,x,y)
        if we is None:
            return
        a,b,x,y = we
        weier = Weierstrass(a,b)
        if point is None:
            return weier
        return Point(weier,x,y)

    


class Edwards(Curve):
    def __init__(self,c: GF,d: GF):
        self.c = c
        self.d = d
        self.field = c.parent()

    def addition(self,point1,point2):
        spoint = edwar_sum(self.c,self.d,point1.x,point1.y,point2.x,point2.y)
        return Point(self,spoint[0],spoint[1])

    def is_infinity(self,x,y,z):
        return x==0 and y==self.c

    def infinity(self):
        return Point(self,self.field(0),self.field(self.c))

    def form(self):
        return "Edwards"

    def negative(self,point):
        return Point(self,-point.x,point.y,point.z)

    def __repr__(self):
        return f"Edwards curve x^2+y^2={self.c}^2(1+{self.d}x^2y^2) over F_{self.field.order()}"

    def check_point(self,x,y,z=1):
        return x**2+y**2==self.c**2*(1+self.d*x**2*y**2)

    def lift_y(self,y):
        x = PolynomialRing(self.field,'x').gen()
        roots = (x**2+y**2-self.c**2*(1+self.d*x**2*y**2)).roots()
        if roots==[]:
            raise NoPoint("No such point")
        return Point(self,roots[0][0],y)

    def __iter__(self):
        try:
            for y in self.field:
                try:
                    point = self.lift_y(y)
                    yield point
                    if not point.x==0: yield -point
                except NoPoint:
                    continue
        except GeneratorExit:
            pass

    def to_weierstrass(self,point=None):
        x,y = (self.field(1),self.field(2)) if point is None else (point.x,point.y)
        we = edwar_to_weier(self.c,self.d,x,y)
        if we is None:
            return
        a,b,x,y = we
        weier = Weierstrass(a,b)
        if point is None:
            return weier
        return Point(weier,x,y)




class Point:
    def __init__(self,curve: Curve,x: GF, y: GF, z: GF=1):
        self.x = x
        self.y = y
        self.z = z
        self.field = x.parent()
        self.curve = curve
        if not self.curve.check_point(x,y,z): raise NoPoint("Point is not on the curve")

    def is_infinity(self):
        return self.curve.is_infinity(self.x,self.y,self.z)

    def __add__(self,other):
        return self.curve.addition(self,other)

    def __rmul__(self,scalar):
        if scalar<0:
            return (-scalar)*(-self)
        accumulator = self.curve.infinity()
        temp = self
        while scalar>0:
            if scalar%2==1:
                accumulator=accumulator+temp
            temp = temp+temp
            scalar>>=1
        return accumulator

    def __eq__(self,other):
        return self.curve==other.curve and self.x==other.x and self.y==other.y and self.z==other.z

    def __repr__(self):
        return f"({self.x},{self.y},{self.z})"

    def __neg__(self):
        return self.curve.negative(self)

    def __sub__(self,other):
        return self+(-other)

    def to_weierstrass(self):
        if self.curve.form()=="Weierstrass":
            return self
        return self.curve.to_weierstrass(self)

    def sage_point(self):
        assert self.curve.form()=="Weierstrass"
        return self.curve.sage_curve(self.x,self.y,self.z)


def test_weierstrass():
    F = GF(101)
    a,b = F(1),F(2)
    W = Weierstrass(a,b)
    I = W.infinity()
    x,y,z = I.x,I.y,I.z
    assert (x,y,z)==(0,1,0)
    assert I.is_infinity()
    P = Point(W,F(77),F(30))
    P6 = 6*P
    assert (P6.x,P6.y,P6.z)==(F(18),F(14),F(1))
    Q = Point(W,F(17),F(36))
    R = P+Q
    assert (R.x,R.y,R.z)==(F(6),F(74),F(1))
    assert (P+I)==P

def test_montgomery():
    F = GF(101)
    a,b = F(49),F(51)
    M = Montgomery(a,b)
    I = M.infinity()
    x,y,z = I.x,I.y,I.z
    assert (x,y,z)==(0,1,0)
    assert I.is_infinity()
    P = Point(M,F(39),F(15))
    P6 = 6*P
    assert (P6.x,P6.y,P6.z)==(F(60),F(7),F(1))
    Q = Point(M,F(9),F(18))
    R = P+Q
    assert (R.x,R.y,R.z)==(F(54),F(37),F(1))
    assert (P+I)==P

def test_twistededwards():
    F = GF(101)
    a,d = F(1),F(94)
    T = TwistedEdwards(a,d)
    I = T.infinity()
    x,y,z = I.x,I.y,I.z
    assert (x,y)==(0,1)
    assert I.is_infinity()
    P = Point(T,F(51),F(21))
    P6 = 6*P
    assert (P6.x,P6.y,P6.z)==(F(39),F(13),F(1))
    Q = Point(T,F(43),F(6))
    R = P+Q
    assert (R.x,R.y,R.z)==(F(97),F(23),F(1))
    assert (P+I)==P

def test_edwards():
    F = GF(101)
    c,d = F(36),F(55)
    T = Edwards(c,d)
    I = T.infinity()
    x,y,z = I.x,I.y,I.z
    assert (x,y)==(0,c)
    assert I.is_infinity()
    P = Point(T,F(29),F(100))
    P6 = 6*P
    assert (P6.x,P6.y,P6.z)==(F(93),F(85),F(1))
    Q = Point(T,F(46),F(17))
    R = P+Q
    assert (R.x,R.y,R.z)==(F(10),F(80),F(1))
    assert (P+I)==P




if __name__=="__main__":
    F = GF(1009)
    w_vec = {'a':F(866),'b':F(208),'P':(F(353),F(449)),'Q':(F(924),F(356)),'R':(F(542),F(665))}
    m_vec = {'a':F(32),'b':F(733),'P':(F(98),F(183)),'Q':(F(915),F(626)),'R':(F(402),F(98))}
    tw_vec = {'a':F(519),'d':F(636),'P':(F(188),F(480)),'Q':(F(593),F(435)),'R':(F(766),F(677))}
    ed_vec = {'c':F(480),'d':F(141),'P':(F(188),F(348)),'Q':(F(593),F(946)),'R':(F(766),F(62))}
    E = EllipticCurve(F,[w_vec['a'],w_vec['b']])
    P = E(*w_vec['P'])
    Q = E(*w_vec['Q'])
    R = E(*w_vec['R'])
    assert P+Q==R

    # test points
    x,y = PolynomialRing(F,['x','y']).gens()
    for point in ['P','Q','R']:
        x,y = w_vec[point]
        assert y**2==x**3+w_vec['a']*x+w_vec['b']
        x,y = m_vec[point]
        assert m_vec['b']*y**2==x**3+m_vec['a']*x**2+x
        x,y = tw_vec[point]
        assert tw_vec['a']*x**2+y**2==1+tw_vec['d']*x**2*y**2
        x,y = ed_vec[point]
        assert x**2+y**2==ed_vec['c']**2*(1+ed_vec['d']*x**2*y**2)

    # test summation formulae
    test_edwar(P,Q,R)
    test_twedw(P,Q,R)
    test_mont(P,Q,R)


    # weier to others
    for point in ['P','Q','R']:
        assert (m_vec['a'],m_vec['b'],*m_vec[point])==weier_to_mont(w_vec['a'],w_vec['b'],*w_vec[point])
        assert (tw_vec['a'],tw_vec['d'],*tw_vec[point])==weier_to_twedw(w_vec['a'],w_vec['b'],*w_vec[point])
        assert (ed_vec['c'],ed_vec['d'],*ed_vec[point])==weier_to_edwar(w_vec['a'],w_vec['b'],*w_vec[point])

    # mont to others
    for point in ['P','Q','R']:
        assert (w_vec['a'],w_vec['b'],*w_vec[point])==mont_to_weier(m_vec['a'],m_vec['b'],*m_vec[point])
        assert (tw_vec['a'],tw_vec['d'],*tw_vec[point])==mont_to_twedw(m_vec['a'],m_vec['b'],*m_vec[point])
        assert (ed_vec['c'],ed_vec['d'],*ed_vec[point])==mont_to_edwar(m_vec['a'],m_vec['b'],*m_vec[point])

    # twedw to others
    for point in ['P','Q','R']:
        assert (w_vec['a'],w_vec['b'],*w_vec[point])==twedw_to_weier(tw_vec['a'],tw_vec['d'],*tw_vec[point])
        assert (m_vec['a'],m_vec['b'],*m_vec[point])==twedw_to_mont(tw_vec['a'],tw_vec['d'],*tw_vec[point])
        assert (ed_vec['c'],ed_vec['d'],*ed_vec[point])==twedw_to_edwar(tw_vec['a'],tw_vec['d'],*tw_vec[point])

    for point in ['P','Q','R']:
        assert (w_vec['a'],w_vec['b'],*w_vec[point])==edwar_to_weier(ed_vec['c'],ed_vec['d'],*ed_vec[point])
        assert (m_vec['a'],m_vec['b'],*m_vec[point])==edwar_to_mont(ed_vec['c'],ed_vec['d'],*ed_vec[point])
        assert (tw_vec['a'],tw_vec['d'],*tw_vec[point])==edwar_to_twedw(ed_vec['c'],ed_vec['d'],*ed_vec[point])

    # test curve classes
    test_weierstrass()
    test_montgomery()
    test_twistededwards()
    test_edwards()
