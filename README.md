Welcome! You have found yourself looking at my personal code repository.

If you're a TA grading my homework, I'm a little bit sorry for you, because this isn't all terribly well-documented. But you probably at least know what functions it is that I used.

## Can I use this for whatever I want?
If you're one of my fellow students (or not), feel free to use this most of this code as you wish. The license, in case you're curious, is GPL-3.0, which means that anyone can copy this and use it as they like, as long as you are willing to put the GPL-3.0 license on the code before you give it to anyone else.

The exceptions to "copy and use how you like" are the following:

Files in the folder `knotts` are not my work; Dr. Thomas Knotts wrote them, and he included licensing information in each of their headers. They use values from the DIPPR thermodynamic database, which you must cite if you publish anything using those values. The classes `muc.wwat`, `muc.wben`, and `muc.wair` use these values.

The code in `steamwrap.py` was not written by myself. 

## How do I use it?

First, know that I make heavy use of a Python module called `pint`. You can find it at pint.readthedocs.io/en/0.16/index.html . It implements units; my code all uses a `pint.UnitRegistry()` object, named `un`. Anytime you see code like `.to("some unit")` or `.magnitude` it is messing around with the units, particularly when they don't play nicely with scipy or numpy. Lots of my functions have an optional argument called either `pint` or `pint_strip`; by default they are set to play nicely with pint units, but you can turn it off by passing `pint=False` or `pint_strip=False`. Unfortunately not all of my code is so gracious, but a lot of it shouldn't depend directly on the units.

My code is also dependent on scipy and numpy, to some extent. I have tried to minimize the outright dependence where possible, but it would be smart of you to make sure you're using the latest versions of numpy and scipy anyway.

The module entitled `muc` (short for "modules, units, constants") is at the core of my work, in the sense that it manages the unit registry for everything else. I always import it. Then, I will use other sections of the code as necessary. Everything does depend on that `muc` module, though, so don't leave it out. It also includes some convenience functions I use a lot, particularly `valprint` (which just prints variables with units nicely) and `list_unit` (which takes a list of values all having the same dimensionality, and converts it to a numpy array with appropriate units).

## What's in here?

`muc` has a number of convenience functions--things I do a lot that I got tired of typing every time. This includes:
* a number of universal constants with `pint` units
* a pretty-printing function called `valprint`
* a similar one called `dictprint` for dictionaries, which wraps `valprint`
* `CtoK`, a function which takes a number, makes it a celsius temperature, and converts to kelvin
* Finite difference derivative, and Explicit Euler integration
 a contains some pint wrappers on Dr. Knott's work, specifically the air, water, and benzene properties files. These wrappers are accessible as `muc.wwat`, `muc.wair`, and `muc.wben`, respectively.

`mixrxn` mainly implements `Mixture` objects, which can represent streams on a physical process. They implement a range of convenient functions, including converting from mass flow rates to mole fractions (or vice versa), a `Convert` method which takes an inlet stream and returns a stream with a given conversion of a given species, and other such helpful things. Courtesy of my newbie self, the code is at least partially documented. `Reaction` objects are designed to play nice with `Mixture` objects, especially for the `Convert` method mentioned.

`thermo` contains a wide range of thermodynamic calculations, the result of a semester in ChEn 373 Thermodynamics. Most of them are object-oriented. Notable examples include EOS calculations (with support for all major cubic EOS and an easy function for calculating pure species fugacity), VLE and other phase equilibria, and reaction equilibrium, complete with activity model functions.

`stats_help` has a function to help me remember how to calculate a confidence interval for the mean of a dataset, but mostly serves to encapsulate some code for calculating confidence and prediction bands for a nonlinear regression.

`spxn` is composed mostly of functions I wrote for homework in ChEn 476, Separations. It also contains some thermodynamic equilibrium data taken from the fourth edition of Wankat's *Separation Process Engineering*, which you should cite if you use the data.


`steamwrap`, by Nathan Barrett, is essentially a `pint` wrapper on jjgomera's IAPWS implementation at github.com/jjgomera/iapws.

## Cool! But it's not working for me.

Shoot me an email and I'll probably respond pretty quickly, or text me if you know me in person and have my phone number. I have no idea how fast I will see issues that get raised here, but you could always try that, too.
