Porthos
=======

***A package of Pari-GP routines for Gaussian quadrature, orthogonal polynomials and related issues, inspired by Walter Gautschi's ORTHPOL.***

Why Pari-GP?
------------

-   Analytical calculations if quantities are exact
    -   rational arithmetic
    -   polynomials and rational functions handled algebraically
-   Numerical calulations when needed or desired
    -   unlimited precision
-   Interactive read-evaluate-print loop
-   Free and libre
-   Well supported
    -   long pedigree
    -   mailing-list
    -   user's manual and tutorials
    -   C subroutines possible
-   Integrated into SAGE

Why another package if Gautschi's is so good?
---------------------------------------------

-   ORTHPOL is not written in Pari-GP
-   Independent languages and software are always useful for verification
-   ORTHPOL and Porthos intersect but neither is a subset of the other

Why the name Porthos?
---------------------

In honour of a fictional 17th-century musketeer from the land of Legendre, Laplace, Laguerre, Hermite, Cauchy, Darboux and Padé, where briefly also Stieltjes lived and died, Henri Cohen conceived Pari-GP, and Bill Allombert and Karim Belabas are currently maintaining it, who liked to just do things first and think about them later.

Quick start
-----------

`porthos54.gp` works with Pari-GP before 2.4.

`porthos60.gp` is almost the same, but works with Pari-GP from 2.4 onwards.

`porthos.gp` is the current version, which at present might work with an earlier release but is developed under Pari 2.7.

Fire up your `gp` with the appropriate package. You will see a message saying:

    <PORTHOS 0.60>  Novice users should type 'help()'.

Do that. You will see:

    help(1): Classified list of functions
    help(2): Types
    help(3): Pseudo-constructors
    help(4): Stieltjes function
    help(5): Recursion coefficients, nodes, gaps and weights
    help(6): Orthogonal expansions and modified moments
    help(7): Left endpoint zero
    help(8): Exact vs numerical calculations
    help(9): Global variables

Porthos does not quite have objects, but it uses Pari-GP's rich variety of recursively structured types to enable you to recognize its five kinds of pseudo-objects by their shape.

### Calculate a 10-point Gaussian quadature rule

    XW(10)

### Find the constant in the leading error term of a quadrature rule

Hardy's rule used to be popular in hand computation days because its coeffients on integer ordinates are exact decimals. When iterated, the pattern of weights is `[0.28, 1.62, 0, 2.2, 0, 1.62, 0.56]`, giving rise to the mnemonic *to remember a number is, to be, a strong or hardy person*. Taking just one panel, the last weight is halved, and on the interval \[-1,1\], everything is divided by 3.

    hardy = [-3, 28/100; -2, 162/100; 0, 220/100; 2, 162/100; 3, 28/100]/3

Method: compute the modified moments with respect to the Legendre polynomials.

    MU(hardy,AB(20))

The first seven terms of the result are

    [2, 0, 0, 0, 0, 0, -4/945]

Thus the degree is 5 (Legendre polynomials 1 to 5 are integrated to 0) and the error constant on the sixth-degree term is -4/945. By way of contrast, let's do the same for the 3-point Gaussian rule. We know this rule also has degree 5, so we go straight for element number 7.

    k6 = MU(XW(3),AB(12))[7]

The answer is -0.045714285714285714285..., clearly a rational number. Pari-GP can work out which it is.

    bestappr(k6)

This is -8/175, almost 11 times as large. A good example of [Brass and Schmeisser's remark](https://www.academia.edu/19806472/Error_estimates_for_interpolatory_quadrature_formulae): "Amongst all quadrature formulae with positive coefficients and fixed order m the Gauss type formulae are worst."
