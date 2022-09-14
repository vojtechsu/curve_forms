---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.8
  kernelspec:
    display_name: SageMath 9.0
    language: sage
    name: sagemath
---

<!-- #region -->
Assume $\text{char}\neq 2,3$


#### (Short) Weierstrass
$y^2 = x^3+ax+b$, 

$4a^3+27b^2\neq 0$

https://en.wikipedia.org/wiki/Elliptic_curve
#### Montgomery

$by^2 = x^3+ax^2+x$,

$b(a^2-4)\neq 0 $

https://en.wikipedia.org/wiki/Montgomery_curve
#### Twisted Edwards

$ax^2+y^2=1+dx^2y^2$,

$(a-d)ad\neq 0 $

https://en.wikipedia.org/wiki/Twisted_Edwards_curve
#### Edwards

$x^2+y^2=c^2(1+dx^2y^2)$,

$cd(1-c^4d)\neq 0 $

https://en.wikipedia.org/wiki/Edwards_curve
#### Hessian

$x^3+y^3+1=3Dxy$,

$D^3\neq 1$

https://en.wikipedia.org/wiki/Hessian_form_of_an_elliptic_curve

#### Legendre

$y^2 = x(x-1)(x-\lambda)$,

$\lambda\neq 0,1$

https://www.math.rug.nl/~top/Legendre.pdf
<!-- #endregion -->

### Transformations

|   	|  Weier 	|   Montg 	|  TwEdw 	|  Edwar  	|  Hessi 	|   Legen	|
|---	|---	|---	|---	|----	|---	|---	|
|   Weier	|   -	| [iso](#weier_to_mont) 	|   [be](#weier_to_twedw)	|  [be](#weier_to_edwar)  	|   	|   	|
|   Montg	|   [iso](#mont_to_weier)	|   - 	|  [be](#mont_to_twedw) 	|   [be](#mont_to_edwar)	|   	|   	|
|   TwEdw	|   [be](#twedw_to_weier)	|  [be](#twedw_to_mont) 	|   -	|  [iso](#twedw_to_edwar)  	|   	|   	|
|   Edwar	|  [be](#edwar_to_weier) 	|  [be](#edwat_to_mont) 	|  [iso](#edwar_to_twedw) 	|   -  	|   	|   	|
|   Hessi	|   	|   	|   	|    	|   -	|   	|
|   Legen	|   	|   	|   	|    	|   	|   - 	|

iso=isomorphism, be=birational equivalence

<!-- #region -->
#### Isomorphism from Weirstrass to Montgomery
<a id='weier_to_mont'></a>
$y^2 = x^3+ax+b \rightarrow Bv^2 = u^3+Au^2+u$, 

$(x,y)\mapsto (u,v)$

Let $\alpha$ be a root of $x^3+ax+b$ and $s=\sqrt{3\alpha^2+a}^{-1}$.

$A= 3\alpha s$, $B=s$ 

$u = s(x-\alpha)$, $v=sy$


Source: https://en.wikipedia.org/wiki/Montgomery_curve

<!-- #endregion -->

```sage
def weier_to_mont(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    for alpha,_ in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            A, B = 3*alpha*s, s
            u,v = s*(x-alpha),s*y
            return A,B,u,v
```

<!-- #region -->
####  Isomorphism from Montgomery to Weirstrass

$by^2 = x^3+ax^2+x \rightarrow v^2 = u^3+Au+B$, 

$(x,y)\mapsto (u,v)$


$A= \frac{3-a^2}{3b^2}$, $B=\frac{2a^3-9a}{27b^3}$ 

$u = \frac{3x+a}{3b}$, $v=\frac{y}{b}$

Source: https://en.wikipedia.org/wiki/Montgomery_curve
<a id='mont_to_weier'></a>
<!-- #endregion -->

```sage
def mont_to_weier(a,b,x,y):
    A, B = (3-a**2)/(3*b**2), (2*a**3-9*a)/(27*b**3)
    u,v = (3*x+a)/(3*b),y/b
    return A,B,u,v
```

#### Birational equivalence from Weierstrass to Twisted Edwards

$y^2 = x^3+ax+b \rightarrow Au^2+v^2=1+Du^2v^2$, 

$(x,y)\mapsto (u,v)$

Let $\alpha$ be a root of $x^3+ax+b$ and $s=\sqrt{3\alpha^2+a}^{-1}$.

$A= \frac{3\alpha s+2}{s}$, $D=\frac{3\alpha s-2}{s}$ 

$u = \frac{x-\alpha}{y}$, $v=\frac{s(x-\alpha)-1}{s(x-\alpha)+1}$

Source: composition of other maps
<a id='weier_to_twedw'></a>

```sage
def weier_to_twedw(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    for alpha,_ in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            A = 3*alpha+2/s
            D = 3*alpha-2/s
            u,v = (x-alpha)/y,(s*(x-alpha)-1)/(s*(x-alpha)+1)
            return A,D,u,v
```

<!-- #region -->
#### Birational equivalence from Twisted Edwards to Weierstrass

$ax^2+y^2=1+dx^2y^2 \rightarrow v^2 = u^3+Au+B$, 

$(x,y)\mapsto (u,v)$

$A= -\frac{a^2+14ad+d^2}{48}$, $B=\frac{(a+d)(-a^2+34ad-d^2)}{864}$ 

$u = \frac{5a+ay-5dy-d}{12-12y}$, $v=\frac{a+ay-dy-d}{4x-4xy}$


The map is a morphism (defined everywhere).

Source: composition of other maps
<a id='twedw_to_weier'></a>
<!-- #endregion -->

```sage
def twedw_to_weier(a,d,x,y):
    A = -(a**2+14*d*a+d**2)/48
    B = (a+d)*(-a**2+34*a*d-d**2)/864
    u = (5*a+a*y-5*d*y-d)/(12-12*y)
    v = (a+a*y-d*y-d)/(4*x-4*x*y)
    return A,B,u,v
```

#### Birational equivalence from Montgomery to Twisted Edwards

$by^2 = x^3+ax^2+x \rightarrow Au^2+v^2=1+Du^2v^2$, 

$(x,y)\mapsto (u,v)$

$A= \frac{a+2}{b}$, $D=\frac{a-2}{b}$ 

$u = \frac{x}{y}$, $v=\frac{x-1}{x+1}$

Source: https://en.wikipedia.org/wiki/Montgomery_curve
<a id='mont_to_twedw'></a>

```sage
def mont_to_twedw(a,b,x,y):
    A,D = (a+2)/b,(a-2)/b
    u,v = x/y,(x-1)/(x+1)
    return A,D,u,v
```

####  Birational equivalence from Twisted Edwards to Montgomery

$ax^2+y^2=1+dx^2y^2 \rightarrow Bv^2 = u^3+Au^2+u$, 

$(x,y)\mapsto (u,v)$

$A= \frac{2(a+d)}{a-d}$, $B=\frac{4}{a-d}$ 

$u = \frac{1+y}{1-y}$, $v=\frac{1+y}{(1-y)x}$

The map is a morphism (defined everywhere).

Source: https://en.wikipedia.org/wiki/Montgomery_curve
<a id='twedw_to_mont'></a>

```sage
def twedw_to_mont(a,d,x,y):
    A,B = 2*(a+d)/(a-d),4/(a-d)
    u,v = (1+y)/(1-y),(1+y)/((1-y)*x)
    return A,B,u,v
```

####  Isomorphism from Twisted Edwards to Edwards (over quadratic extension)

$ax^2+y^2=1+dx^2y^2 \rightarrow u^2+v^2=C^2(1+Du^2v^2)$, 

$(x,y)\mapsto (u,v)$

$C= \frac{1}{\sqrt{a}}$, $D=ad$ 

$u = x$, $v=\frac{y}{\sqrt{a}}$

Source: I thought really hard
<a id='twedw_to_edwar'></a>

```sage
def twedw_to_edwar(a,d,x,y):
    p = a.parent().order()
    sa = sqrt(1/a) if ZZ(sqrt(1/a))<ZZ(p-sqrt(1/a)) else -sqrt(1/a) # take the smallest square root
    C,D = sa,d*a
    u,v = x,y*C
    return C,D,u,v
```

####  Isomorphism from Edwards to Twisted Edwards

$x^2+y^2=c^2(1+dx^2y^2) \rightarrow Au^2+v^2=1+Du^2v^2$, 

$(x,y)\mapsto (u,v)$

$A= \frac{1}{c^2}$, $D=c^2d$ 

$u = x$, $v=\frac{y}{c}$

Source: I thought really hard
<a id='edwar_to_twedw'></a>

```sage
def edwar_to_twedw(c,d,x,y):
    u,v = x,y/c
    A,D = 1/c**2, c**2*d
    return A,D,u,v
```

####  Birational equivalence from Edwards to Weierstrass

$x^2+y^2=c^2(1+dx^2y^2) \rightarrow v^2=u^3+Au^2+B$, 

$(x,y)\mapsto (u,v)$

$A= -\frac{1+14dc^4+c^8d^2}{48c^4}$, $B=\frac{(1+c^4d)(-1+34dc^4-c^8d^2)}{864c^6}$ 

$u = \frac{5c+y-5c^4dy-c^5d}{12c^3-12yc^2}$, $v=\frac{c+y-5dyc^4-c^5d}{4xc^3-4xyc^2}$

Source: composition of other maps
<a id='edwar_to_weier'></a>

```sage
def edwar_to_weier(c,d,x,y):
    u,v = (5*c+y-5*c**4*d*y-c**5*d)/(12*c**3-12*y*c**2),(c+y-d*y*c**4-c**5*d)/(4*x*c**3-4*x*y*c**2)
    A,B = -(1+14*d*c**4+c**8*d**2)/(48*c**4), (1+c**4*d)*(-1+34*d*c**4-c**8*d**2)/(864*c**6)
    return A,B,u,v
```

#### Birational equivalence from Weierstrass to Edwards

$y^2 = x^3+ax+b \rightarrow u^2+v^2=C^2(1+Du^2v^2)$, 

$(x,y)\mapsto (u,v)$

Let $\alpha$ be a root of $x^3+ax+b$ and $s=\sqrt{3\alpha^2+a}^{-1}$.

$A= \sqrt{\frac{s}{3s\alpha+2}}$, $D=-4a-3\alpha^2$ 

$u = \frac{x-\alpha}{y}$, $v=\frac{s(x-\alpha)-1}{s(x-\alpha)+1}\cdot\sqrt{\frac{s}{3s\alpha+2}}$

Source: composition of other maps
<a id='weier_to_edwar'></a>

```sage
def weier_to_edwar(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    p = a.parent().order()
    for alpha,ex in (z**3+a*z+b).roots():
        if (3*alpha**2+a).is_square():
            s = 1/(sqrt(3*alpha**2+a))
            t = sqrt(s/(3*s*alpha+2))
            t = t if ZZ(t)<ZZ(p-t) else p-t # take the smallest square root
            C, D = t,-4*a-3*alpha**2
            u,v = (x-alpha)/y,(s*(x-alpha)-1)/(s*(x-alpha)+1)*t
            return C,D,u,v
```

#### Birational equivalence from Montgomery to Edwards

$by^2 = x^3+ax^2+x \rightarrow u^2+v^2=C^2(1+Du^2v^2)$, 

$(x,y)\mapsto (u,v)$

$C= \sqrt{\frac{b}{a+2}}$, $D=\frac{a^2-4}{b^2}$ 

$u = \frac{x}{y}$, $v=\frac{x-1}{x+1} \sqrt{\frac{b}{a+2}}$

Source: composition of other maps
<a id='mont_to_edwar'></a>

```sage
def mont_to_edwar(a,b,x,y):
    C,D = sqrt(b/(a+2)),(a**2-4)/(b**2)
    u,v = x/y,(x-1)/(x+1)*sqrt(b/(a+2))
    return C,D,u,v
```

####  Birational equivalence from Edwards to Montgomery

$x^2+y^2=c^2(1+dx^2y^2) \rightarrow Bv^2 = u^3+Au^2+u$, 

$(x,y)\mapsto (u,v)$

$A= \frac{2+2c^4d}{1-c^4d}$, $B=\frac{4c^2}{1-c^4d}$ 

$u = \frac{c+y}{c-y}$, $v=\frac{c+y}{(c-y)x}$

Source: composition of other maps
<a id='edwar_to_mont'></a>

```sage
def edwar_to_mont(c,d,x,y):
    A,B = (2+2*c**4*d)/(1-c**4*d),(4*c**2)/(1-c**4*d)
    u,v = (c+y)/(c-y),(c+y)/((c-y)*x)
    return A,B,u,v
```

```sage

```

```sage

```
