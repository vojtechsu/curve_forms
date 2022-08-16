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


#### Weierstrass
$y^2 = x^3+ax+b$, 

$4a^3+27b^2\neq 0$

https://en.wikipedia.org/wiki/Elliptic_curve
#### Montgomery

$by^2 = x^3+ax^2+x$,

$b(a^2-4)\neq =0 $

https://en.wikipedia.org/wiki/Montgomery_curve
#### Twisted Edwards

$ax^2+y^2=1+dx^2y^2$

https://en.wikipedia.org/wiki/Twisted_Edwards_curve
#### Edwards

$x^2+y^2=1+dx^2y^2$

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
|   Weier	|   -	| [ref](#weier_to_mont) 	|   [ref](#weier_to_twedw)	|    	|   	|   	|
|   Montg	|   [ref](#mont_to_weier)	|   - 	|  [ref](#mont_to_twedw) 	|   	|   	|   	|
|   TwEdw	|   [ref](#twedw_to_weier)	|  [ref](#twedw_to_mont) 	|   -	|    	|   	|   	|
|   Edwar	|   	|   	|   	|   -  	|   	|   	|
|   Hessi	|   	|   	|   	|    	|   -	|   	|
|   Legen	|   	|   	|   	|    	|   	|   - 	|


#### Weirstrass to Montgomery
<a id='weier_to_mont'></a>
$y^2 = x^3+ax+b \rightarrow Bv^2 = u^3+Au^2+u$, 

$(x,y)\mapsto (u,v)$

Let $\alpha$ be a root of $x^3+ax+b$ and $s=\sqrt{3\alpha^2+a}^{-1}$.

$A= 3\alpha s$, $B=s$ 

$u = s(x-\alpha)$, $v=sy$

https://en.wikipedia.org/wiki/Montgomery_curve


```sage
def weier_to_mont(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    try:
        alpha = (z**3+a*z+b).roots()[0][0]
        s = 1/(sqrt(3*alpha**2+a))
    except:
        "No Montgomery form"
    A, B = 3*alpha*s, s
    u,v = s*(x-alpha),s*y
    return A,B,u,v
```

<!-- #region -->
####  Montgomery to Weirstrass

$by^2 = x^3+ax^2+x \rightarrow v^2 = u^3+Au+B$, 

$(x,y)\mapsto (u,v)$


$A= \frac{3-a^2}{3b^2}$, $B=\frac{2a^3-9a}{27b^3}$ 

$u = \frac{3x+a}{b}$, $v=\frac{y}{b}$

https://en.wikipedia.org/wiki/Montgomery_curve
<a id='mont_to_weier'></a>
<!-- #endregion -->

```sage
def mont_to_weier(a,b,x,y):
    A, B = (3-a**2)/(3*b**2), (2*a**3-9*a)/(27*b**2)
    u,v = (3*x+a)/b,y/b
    return A,B,u,v
```

#### Weierstrass to Twisted Edwards

$y^2 = x^3+ax+b \rightarrow Au^2+v^2=1+Du^2v^2$, 

$(x,y)\mapsto (u,v)$

Let $\alpha$ be a root of $x^3+ax+b$ and $s=\sqrt{3\alpha^2+a}^{-1}$.

$A= \frac{3\alpha s+2}{s}$, $D=\frac{3\alpha s-2}{s}$ 

$u = \frac{x-\alpha}{y}$, $v=\frac{s(x-\alpha)-1}{s(x-\alpha)+1}$

<a id='weier_to_twedw'></a>

```sage
def weier_to_twedw(a,b,x,y):
    z = PolynomialRing(a.parent(),'z').gen()
    try:
        alpha = (z**3+a*z+b).roots()[0][0]
        s = 1/(sqrt(3*alpha**2+a))
    except:
        "No Twisted Edwards form"
    A = 3*alpha+2/s
    D = 3*alpha-2/s
    u,v = (x-alpha)/y,(s*(x-alpha)-1)/(s*(x-alpha)+1)
    return A,D,u,v
```

#### Twisted Edwards to Weierstrass

$au^2+v^2=1+du^2v^2 \rightarrow v^2 = u^3+Au+B$, 

$(x,y)\mapsto (u,v)$

$A= -\frac{a^2+14ad+d^2}{48}$, $B=\frac{(a+d)(-a^2+34ad-d^2)}{864}$ 

$u = \frac{5a+ay-5dy-d}{12-12y}$, $v=\frac{a+ay-dy-d}{4x-4xy}$

<a id='twedw_to_weier'></a>

```sage
def twedw_to_weier(a,d,x,y):
    A = -(a**2+14*a*d+d**2)/48
    B = (a+d)*(-a**2+34*a*d-d**2)/864
    u = (5*a+a*y-5*d*y-d)/(12-12*y)
    v = (a+a*y-d*y-d)/(4*x-4*x*y)
    return A,B,u,v
```

#### Montgomery to Twisted Edwards

$by^2 = x^3+ax^2+x \rightarrow Au^2+v^2=1+Du^2v^2$, 

$(x,y)\mapsto (u,v)$

$A= \frac{a+2}{b}$, $D=\frac{a-2}{b}$ 

$u = \frac{x}{y}$, $v=\frac{x-1}{x+1}$

<a id='mont_to_twedw'></a>

```sage
def mont_to_twedw(a,b,x,y):
    A,B = (a+2)/b,(a-2)/b
    u,v = x/y,(x-1)/(x+1)
    return A,B,u,v
```

####  Twisted Edwards to Montgomery

$au^2+v^2=1+du^2v^2 \rightarrow Bv^2 = u^3+Au^2+u$, 

$(x,y)\mapsto (u,v)$

$A= \frac{2(a+d)}{a-d}$, $B=\frac{4}{a-d}$ 

$u = \frac{1+y}{1-y}$, $v=\frac{1+y}{(1-y)x}$

<a id='twedw_to_mont'></a>

```sage
def twedw_to_mont(a,b,x,y):
    A,B = 2*(a+d)/(a-d),4/(a-d)
    u,v = (1+y)/(1-y),(1+y)/((1-y)*x)
    return A,B,u,v
```

```sage

```

```sage

```
